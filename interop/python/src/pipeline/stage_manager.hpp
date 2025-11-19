#ifndef LAHUTA_BINDINGS_STAGE_MANAGER_HPP
#define LAHUTA_BINDINGS_STAGE_MANAGER_HPP

#include <chrono>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include "analysis/contacts/computation.hpp"
#include "analysis/contacts/provider.hpp"
#include "analysis/system/model_pack_task.hpp"
#include "interactions.hpp"
#include "pipeline/compute/parameters.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/dynamic/run_metrics.hpp"
#include "pipeline/process_pool.hpp"
#include "pipeline/process_task.hpp"
#include "pipeline/thread_task.hpp"

// clang-format off
namespace py = pybind11;
namespace lahuta::bindings {
using namespace lahuta::sources;
using namespace lahuta::pipeline::dynamic;

namespace {

inline InteractionTypeSet normalize_interaction_argument(py::handle obj) {
  if (!obj || obj.is_none()) return InteractionTypeSet::all();
  InteractionTypeSet set;
  add_python_interactions(set, obj);
  if (set.empty()) {
    throw py::value_error("Interaction type selection cannot be empty");
  }
  return set;
}

using DefaultExecutor = StageExecutor<StageRunMetrics>;
using NullExecutor    = StageExecutor<NullStageRunMetrics>;

inline void python_stage_executor_tls_cleanup(DefaultExecutor::ThreadLocalState& state) {
  // During interpreter shutdown we cannot safely acquire the GIL. Fall back to direct reset.
  if (!Py_IsInitialized() || Py_IsFinalizing()) {
    DefaultExecutor::clear_thread_local_state(state);
    return;
  }
  try {
    py::gil_scoped_acquire gil;
    DefaultExecutor::clear_thread_local_state(state);
  } catch (const py::error_already_set&) {
    DefaultExecutor::clear_thread_local_state(state);
  } catch (...) {
    DefaultExecutor::clear_thread_local_state(state);
  }
}

inline void python_stage_executor_tls_cleanup_null(NullExecutor::ThreadLocalState& state) {
  if (!Py_IsInitialized() || Py_IsFinalizing()) {
    NullExecutor::clear_thread_local_state(state);
    return;
  }
  try {
    py::gil_scoped_acquire gil;
    NullExecutor::clear_thread_local_state(state);
  } catch (const py::error_already_set&) {
    NullExecutor::clear_thread_local_state(state);
  } catch (...) {
    NullExecutor::clear_thread_local_state(state);
  }
}

inline std::chrono::milliseconds seconds_to_milliseconds(double seconds) {
  using MilliRep = std::chrono::milliseconds::rep;
  if (!std::isfinite(seconds)) {
    throw std::invalid_argument("flush_timeout must be a finite value");
  }
  if (seconds < 0.0) {
    throw std::invalid_argument("flush_timeout must be non-negative");
  }
  const double max_seconds = static_cast<double>(std::numeric_limits<MilliRep>::max()) / 1000.0;
  if (seconds > max_seconds) {
    throw std::invalid_argument("flush_timeout is too large");
  }
  auto duration = std::chrono::duration<double>(seconds);
  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(duration);
  return ms;
}

inline double milliseconds_to_seconds(std::chrono::milliseconds ms) {
  return std::chrono::duration<double>(ms).count();
}

} // namespace

inline void bind_stage_manager(py::module_ &md) {
  DefaultExecutor::set_tls_cleanup_hook(&python_stage_executor_tls_cleanup);
  NullExecutor::set_tls_cleanup_hook(&python_stage_executor_tls_cleanup_null);

  py::enum_<StageManager::ReportingLevel>(md, "ReportingLevel")
    .value("OFF",   StageManager::ReportingLevel::Off)
    .value("BASIC", StageManager::ReportingLevel::Basic)
    .value("DEBUG", StageManager::ReportingLevel::Debug);

  py::class_<ITask, std::shared_ptr<ITask>> task(md, "Task");

  // Bind reusable C++ task that serializes models to ModelRecord payloads
  py::class_<analysis::system::ModelPackTask, ITask, std::shared_ptr<analysis::system::ModelPackTask>>(md, "ModelPackTask")
    .def(py::init<std::string>(), py::arg("channel") = std::string("db"));

  py::class_<StageManager, std::shared_ptr<StageManager>>(md, "StageManager")
    .def(py::init([](std::shared_ptr<sources::IDescriptor> source) {
          return std::make_shared<StageManager>(source);
        }),
        py::arg("source"))
    // Generic task registration
    .def("add_task", [](StageManager& mgr, const std::string& name,
                         const std::vector<std::string>& deps,
                         std::shared_ptr<ITask> task,
                         bool thread_safe) {
          mgr.add_task(name, deps, std::move(task), thread_safe);
        },
        py::arg("name"), py::arg("depends") = std::vector<std::string>{},
        py::arg("task"), py::arg("thread_safe") = true)

    .def("add_contacts", [](StageManager &mgr,
                            const std::string &name,
                            analysis::contacts::ContactProvider provider,
                            py::object interaction_obj,
                            std::optional<std::string> channel,
                            const std::string &out_fmt,
                            bool thread_safe) {
          std::string ch = channel.has_value() ? *channel : name;

          pipeline::compute::ContactsOutputFormat format;
          if      (out_fmt == "json")   format = pipeline::compute::ContactsOutputFormat::Json;
          else if (out_fmt == "text")   format = pipeline::compute::ContactsOutputFormat::Text;
          else if (out_fmt == "binary") format = pipeline::compute::ContactsOutputFormat::Binary;
          else throw std::invalid_argument("out_fmt must be 'json', 'text', or 'binary'");

          InteractionTypeSet interaction_types = normalize_interaction_argument(interaction_obj);

          pipeline::compute::ContactsParams p{};
          p.provider = provider;
          p.type     = interaction_types;
          p.channel  = ch;
          p.format   = format;

          std::vector<std::string> deps = {"topology"};
          mgr.add_computation(name, deps, [label = name, p]() {
            return std::make_unique<analysis::contacts::ContactsComputation>(label, p);
          }, thread_safe);
        },
        py::arg("name"),
        py::arg("provider"), py::arg("interaction_types") = py::none(),
        py::arg("channel") = std::optional<std::string>{},
        py::arg("out_fmt") = std::string("json"),
        py::arg("thread_safe") = true)

    // python tasks
    .def("add_python", [](StageManager &mgr,
                          const std::string &name,
                          const std::vector<std::string> &deps,
                          py::object fn,
                          std::optional<std::string> channel,
                          bool serialize,
                          bool store) {
          // Always pass PipelineContext to the callable. Serialize return values as text or JSON (see PyCallableTask)
          std::optional<std::string> ch = channel.has_value() ? channel : std::optional<std::string>{name};

          auto t = std::make_shared<PyCallableTask>(std::move(fn), name, ch, store, serialize);
          // Preserve item-level parallelism: mark stage as thread-safe.
          mgr.add_task(name, deps, std::move(t), /*thread_safe=*/true);
        },
        py::arg("name"),
        py::arg("depends") = std::vector<std::string>{},
        py::arg("fn"),
        py::arg("channel") = std::optional<std::string>{},
        py::arg("serialize") = true,
        py::arg("store") = true)
    .def("add_python_process", [](StageManager &mgr,
                                  const std::string &name,
                                  const std::vector<std::string> &deps,
                                  const std::string &module,
                                  const std::string &qualname,
                                  py::object serialized_callable,
                                  std::optional<std::string> channel,
                                  bool store) {
          std::optional<std::string> ch = channel.has_value() ? channel : std::optional<std::string>{name};
          auto& sys_params = mgr.get_system_params();
          auto& topo_params = mgr.get_topology_params();
          auto t = std::make_shared<PyProcessTask>(name,
                                                   module,
                                                   qualname,
                                                   std::move(serialized_callable),
                                                   ch,
                                                   store,
                                                   &sys_params,
                                                   &topo_params);
          mgr.add_task(name, deps, std::move(t), /*thread_safe=*/true);
        },
        py::arg("name"),
        py::arg("depends") = std::vector<std::string>{},
        py::arg("module"),
        py::arg("qualname"),
        py::arg("serialized_callable") = py::none(),
        py::arg("channel") = std::optional<std::string>{},
        py::arg("store") = true)

    // sinks
    .def("connect_sink", [](StageManager &mgr,
                             const std::string &channel,
                             std::shared_ptr<IDynamicSink> sink,
                             std::optional<BackpressureConfig> cfg) {
          mgr.connect_sink(channel, std::move(sink), std::move(cfg));
        },
        py::arg("channel"), py::arg("sink"), py::arg("backpressure") = std::nullopt)
    .def("sorted_tasks", &StageManager::sorted_tasks)
    .def("compile",      &StageManager::compile)
    .def("stats", [](StageManager &mgr) {
          py::list out;
          auto v = mgr.stats();
          for (const auto &s : v) {
            py::dict d;
            d["sink_name"]      = s.sink_name;
            d["enqueued_msgs"]  = s.enqueued_msgs;
            d["enqueued_bytes"] = s.enqueued_bytes;
            d["written_msgs"]   = s.written_msgs;
            d["written_bytes"]  = s.written_bytes;
            d["stalled_ns"]     = s.stalled_ns;
            d["drops"]          = s.drops;
            d["queue_msgs"]     = s.queue_msgs;
            d["queue_bytes"]    = s.queue_bytes;
            d["writer_threads"] = s.writer_threads;
            out.append(std::move(d));
          }
          return out;
        })
    .def("run", [](StageManager &mgr, int threads, std::optional<double> flush_timeout) {
          std::optional<std::chrono::milliseconds> previous_timeout;
          if (flush_timeout.has_value()) {
            const auto override_timeout = seconds_to_milliseconds(*flush_timeout);
            previous_timeout = mgr.get_flush_timeout();
            mgr.set_flush_timeout(override_timeout);
          }
          struct FlushTimeoutGuard {
            StageManager& mgr;
            std::optional<std::chrono::milliseconds> previous;
            ~FlushTimeoutGuard() {
              if (previous.has_value()) {
                mgr.set_flush_timeout(*previous);
              }
            }
          } guard{mgr, previous_timeout};
          {
            py::gil_scoped_release release;
            mgr.run(threads);
          }
          py::gil_scoped_acquire acquire;
        },
        py::arg("threads") = 4, py::arg("flush_timeout") = py::none())
    .def("last_run_report", [](StageManager &mgr) -> py::object {
          const auto& opt = mgr.last_report();
          if (!opt.has_value()) return py::none();
          const auto& report = *opt;
          py::dict d;
          d["total_seconds"]     = report.total_seconds;
          d["cpu_seconds"]       = report.cpu_seconds;
          d["io_seconds"]        = report.io_seconds;
          d["ingest_seconds"]    = report.ingest_seconds;
          d["prepare_seconds"]   = report.prepare_seconds;
          d["flush_seconds"]     = report.flush_seconds;
          d["setup_seconds"]     = report.setup_seconds;
          d["compute_seconds"]   = report.compute_seconds;
          d["items_total"]       = report.items_total;
          d["items_processed"]   = report.items_processed;
          d["items_skipped"]     = report.items_skipped;
          d["stage_count"]       = report.stage_count;
          d["threads_requested"] = report.threads_requested;
          d["threads_used"]      = report.threads_used;
          d["all_thread_safe"]   = report.all_thread_safe;
          d["run_token"]         = report.run_token;
          d["metrics_enabled"]   = report.metrics_enabled;
          d["peak_inflight_items"]   = report.peak_inflight_items;
          d["average_queue_depth"]   = report.average_queue_depth;
          d["permit_wait_events"]    = report.permit_wait_events;
          d["permit_wait_total_seconds"] = report.permit_wait_total_seconds;
          d["permit_wait_min_seconds"]   = report.permit_wait_min_seconds;
          d["permit_wait_max_seconds"]   = report.permit_wait_max_seconds;
          d["permit_wait_avg_seconds"]   = report.permit_wait_avg_seconds;
          py::list stage_breakdown;
          for (const auto& st : report.stage_breakdown) {
            py::dict stage_entry;
            stage_entry["label"] = st.label;
            stage_entry["setup_seconds"] = st.setup_seconds;
            stage_entry["compute_seconds"] = st.compute_seconds;
            stage_breakdown.append(std::move(stage_entry));
          }
          d["stage_breakdown"] = std::move(stage_breakdown);
          d["mux_sink_count"]        = report.mux_sink_count;
          d["mux_enqueued_msgs"]     = report.mux_enqueued_msgs;
          d["mux_enqueued_bytes"]    = report.mux_enqueued_bytes;
          d["mux_written_msgs"]      = report.mux_written_msgs;
          d["mux_written_bytes"]     = report.mux_written_bytes;
          d["mux_stall_ns"]          = report.mux_stall_ns;
          d["mux_drops"]             = report.mux_drops;
          d["mux_queue_depth_peak"]  = report.mux_queue_depth_peak;
          d["mux_queue_bytes_peak"]  = report.mux_queue_bytes_peak;
          d["mux_active_writers_total"] = report.mux_active_writers_total;
          d["mux_active_writers_peak"]  = report.mux_active_writers_peak;
          return d;
        })

    // Process pool management
    .def("configure_python_process_pool", [](StageManager&, std::size_t processes) {
          PyProcessPool::instance().configure(processes);
        }, py::arg("processes"))
    .def("shutdown_python_process_pool", [](StageManager&) {
          PyProcessPool::instance().shutdown();
          PyProcessTask::set_concurrency_limit(0);
        })
    .def("set_python_process_concurrency", [](StageManager&, std::size_t limit) {
          PyProcessTask::set_concurrency_limit(limit);
        }, py::arg("limit"))
    .def("set_python_process_timeout", [](StageManager&, double seconds) {
          PyProcessTask::set_timeout(seconds);
        }, py::arg("seconds"))

    .def("set_auto_builtins", [](StageManager& mgr, bool on) { mgr.set_auto_builtins(on); })
    .def("get_auto_builtins", [](StageManager& mgr) { return mgr.get_auto_builtins(); })
    .def("set_reporting_level", [](StageManager& mgr, StageManager::ReportingLevel level) {
          mgr.set_reporting_level(level);
        })
    .def("get_reporting_level", [](StageManager& mgr) {
          return mgr.get_reporting_level();
        })

    .def("set_flush_timeout", [](StageManager& mgr, double timeout_seconds) {
          mgr.set_flush_timeout(seconds_to_milliseconds(timeout_seconds));
        }, py::arg("timeout"))
    .def("get_flush_timeout", [](StageManager& mgr) {
          return milliseconds_to_seconds(mgr.get_flush_timeout());
        })

    // Parameter access for builtins
    .def("get_system_params", [](StageManager &mgr) -> py::dict {
          const auto& params = mgr.get_system_params();
          py::dict result;
          result["is_model"] = params.is_model;
          return result;
        })
    .def("set_system_params", [](StageManager &mgr, py::dict d) {
          auto& params = mgr.get_system_params();
          if (d.contains("is_model")) params.is_model = d["is_model"].cast<bool>();
          mgr.invalidate_compilation();
        })
    .def("get_topology_params", [](StageManager &mgr) -> py::dict {
          const auto& params = mgr.get_topology_params();
          py::dict result;
          result["flags"] = static_cast<int>(params.flags);
          result["atom_typing_method"] = static_cast<int>(params.atom_typing_method);
          return result;
        })
    .def("set_topology_params", [](StageManager &mgr, py::dict d) {
          auto& params = mgr.get_topology_params();
          if (d.contains("flags")) {
            // Handle both int and enum types
            if (py::isinstance<py::int_>(d["flags"])) {
              params.flags = static_cast<TopologyComputation>(d["flags"].cast<int>());
            } else {
              params.flags = d["flags"].cast<TopologyComputation>();
            }
          }
          if (d.contains("atom_typing_method")) {
            if (py::isinstance<py::int_>(d["atom_typing_method"])) {
              params.atom_typing_method = static_cast<AtomTypingMethod>(d["atom_typing_method"].cast<int>());
            } else {
              params.atom_typing_method = d["atom_typing_method"].cast<AtomTypingMethod>();
            }
          }
          mgr.invalidate_compilation();
        })

    // Graph introspection
    .def("describe_graph", &StageManager::describe_graph);
}

} // namespace lahuta::bindings

#endif // LAHUTA_BINDINGS_STAGE_MANAGER_HPP
