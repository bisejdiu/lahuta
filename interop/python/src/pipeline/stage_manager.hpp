/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, '\0'); char* ptr = s.data();
 *   for (char c : std::string_view{"besian"}) *ptr++ = c;
 *   for (char c : std::string_view{"sejdiu"}) *ptr++ = c;
 *   *ptr++ = '@';
 *   for (char c : std::string_view{"gmail.com"}) *ptr++ = c;
 *   return s;
 * }();
 *
 */

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
#include "pipeline/io/backpressure.hpp"
#include "pipeline/metrics/run_metrics.hpp"
#include "pipeline/runtime/manager.hpp"
#include "pipeline/task/compute/parameters.hpp"
#include "pipeline/thread_task.hpp"

namespace py = pybind11;
namespace lahuta::bindings {
namespace P = lahuta::pipeline;

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

using DefaultExecutor = P::StageExecutor<P::StageRunMetrics>;
using NullExecutor    = P::StageExecutor<P::NullStageRunMetrics>;

inline void python_stage_executor_tls_cleanup(DefaultExecutor::ThreadLocalState &state) {
  // During interpreter shutdown we cannot safely acquire the GIL. Fall back to direct reset.
  if (!Py_IsInitialized() || Py_IsFinalizing()) {
    DefaultExecutor::clear_thread_local_state(state);
    return;
  }
  try {
    py::gil_scoped_acquire gil;
    DefaultExecutor::clear_thread_local_state(state);
  } catch (const py::error_already_set &) {
    DefaultExecutor::clear_thread_local_state(state);
  } catch (...) {
    DefaultExecutor::clear_thread_local_state(state);
  }
}

inline void python_stage_executor_tls_cleanup_null(NullExecutor::ThreadLocalState &state) {
  if (!Py_IsInitialized() || Py_IsFinalizing()) {
    NullExecutor::clear_thread_local_state(state);
    return;
  }
  try {
    py::gil_scoped_acquire gil;
    NullExecutor::clear_thread_local_state(state);
  } catch (const py::error_already_set &) {
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
  auto ms       = std::chrono::duration_cast<std::chrono::milliseconds>(duration);
  return ms;
}

inline double milliseconds_to_seconds(std::chrono::milliseconds ms) {
  return std::chrono::duration<double>(ms).count();
}

} // namespace

inline std::vector<P::DataField> bitset_to_fields(P::DataFieldSet set) {
  std::vector<P::DataField> fields;
  auto try_add = [&](P::DataField field) {
    if (set.contains(field)) fields.push_back(field);
  };
  try_add(P::DataField::Metadata);
  try_add(P::DataField::Sequence);
  try_add(P::DataField::Positions);
  try_add(P::DataField::Plddt);
  try_add(P::DataField::Dssp);
  try_add(P::DataField::SequenceView);
  try_add(P::DataField::PositionsView);
  try_add(P::DataField::PlddtView);
  try_add(P::DataField::DsspView);
  return fields;
}

inline P::DataFieldSet fields_to_bitset(const std::vector<P::DataField> &fields) {
  P::DataFieldSet set = P::DataFieldSet::none();
  for (auto field : fields) {
    set |= field;
  }
  return set;
}

inline void bind_stage_manager(py::module_ &md) {
  DefaultExecutor::set_tls_cleanup_hook(&python_stage_executor_tls_cleanup);
  NullExecutor::set_tls_cleanup_hook(&python_stage_executor_tls_cleanup_null);

  py::enum_<P::DataField>(md, "DataField")
      .value("Metadata", P::DataField::Metadata)
      .value("Sequence", P::DataField::Sequence)
      .value("Positions", P::DataField::Positions)
      .value("Plddt", P::DataField::Plddt)
      .value("Dssp", P::DataField::Dssp)
      .value("SequenceView", P::DataField::SequenceView)
      .value("PositionsView", P::DataField::PositionsView)
      .value("PlddtView", P::DataField::PlddtView)
      .value("DsspView", P::DataField::DsspView);

  py::enum_<P::StageManager::ReportingLevel>(md, "ReportingLevel")
      .value("OFF", P::StageManager::ReportingLevel::Off)
      .value("BASIC", P::StageManager::ReportingLevel::Basic)
      .value("DEBUG", P::StageManager::ReportingLevel::Debug);

  py::class_<P::ITask, std::shared_ptr<P::ITask>> task(md, "Task");

  // Bind reusable C++ task that serializes models to ModelRecord payloads
  py::class_<analysis::ModelPackTask, P::ITask, std::shared_ptr<analysis::ModelPackTask>>(md, "ModelPackTask")
      .def(py::init<std::string>(), py::arg("channel") = std::string("db"));

  py::class_<P::StageManager, std::shared_ptr<P::StageManager>>(md, "StageManager")
      .def(py::init([](std::shared_ptr<P::IDescriptor> source) {
             return std::make_shared<P::StageManager>(source);
           }),
           py::arg("source"))
      // Generic task registration
      .def(
          "add_task",
          [](P::StageManager &mgr,
             const std::string &name,
             const std::vector<std::string> &deps,
             std::shared_ptr<P::ITask>
                 task,
             bool thread_safe) { mgr.add_task(name, deps, std::move(task), thread_safe); },
          py::arg("name"),
          py::arg("depends") = std::vector<std::string>{},
          py::arg("task"),
          py::arg("thread_safe") = true)

      .def(
          "add_contacts",
          [](P::StageManager &mgr,
             const std::string &name,
             analysis::ContactProvider provider,
             py::object interaction_obj,
             std::optional<std::string>
                 channel,
             const std::string &out_fmt,
             bool thread_safe) {
            std::string ch = channel.has_value() ? *channel : name;

            P::ContactsOutputFormat format;
            if (out_fmt == "json")
              format = P::ContactsOutputFormat::Json;
            else if (out_fmt == "binary")
              format = P::ContactsOutputFormat::Binary;
            else
              throw std::invalid_argument("out_fmt must be 'json' or 'binary'");

            InteractionTypeSet interaction_types = normalize_interaction_argument(interaction_obj);

            P::ContactsParams p{};
            p.provider = provider;
            p.type     = interaction_types;
            p.channel  = ch;
            p.format   = format;

            std::vector<std::string> deps = {"topology"};
            mgr.add_computation(
                name,
                deps,
                [label = name, p]() { return std::make_unique<analysis::ContactsComputation>(label, p); },
                thread_safe);
          },
          py::arg("name"),
          py::arg("provider"),
          py::arg("interaction_types") = py::none(),
          py::arg("channel")           = std::optional<std::string>{},
          py::arg("out_fmt")           = std::string("json"),
          py::arg("thread_safe")       = true)

      .def(
          "set_task_data_fields",
          [](P::StageManager &mgr, const std::string &name, const std::vector<P::DataField> &fields) {
            mgr.set_task_data_requirements(name, fields_to_bitset(fields));
          },
          py::arg("name"),
          py::arg("fields"))

      .def(
          "get_task_data_fields",
          [](P::StageManager &mgr, const std::string &name) {
            return bitset_to_fields(mgr.get_task_data_requirements(name));
          },
          py::arg("name"))

      // python tasks
      .def(
          "add_python",
          [](P::StageManager &mgr,
             const std::string &name,
             const std::vector<std::string> &deps,
             py::object fn,
             std::optional<std::string>
                 channel,
             bool serialize,
             bool store) {
            // Always pass PipelineContext to the callable. Serialize return values as text or JSON (see
            // PyCallableTask)
            std::optional<std::string> ch = channel.has_value() ? channel : std::optional<std::string>{name};

            auto t = std::make_shared<PyCallableTask>(std::move(fn), name, ch, store, serialize);
            // Preserve item-level parallelism: mark stage as thread-safe.
            mgr.add_task(name, deps, std::move(t), /*thread_safe=*/true);
          },
          py::arg("name"),
          py::arg("depends") = std::vector<std::string>{},
          py::arg("fn"),
          py::arg("channel")   = std::optional<std::string>{},
          py::arg("serialize") = true,
          py::arg("store")     = true)

      // sinks
      .def(
          "connect_sink",
          [](P::StageManager &mgr,
             const std::string &channel,
             std::shared_ptr<P::IDynamicSink>
                 sink,
             std::optional<P::BackpressureConfig>
                 cfg) { mgr.connect_sink(channel, std::move(sink), std::move(cfg)); },
          py::arg("channel"),
          py::arg("sink"),
          py::arg("backpressure") = std::nullopt)
      .def("sorted_tasks", &P::StageManager::sorted_tasks)
      .def("compile", &P::StageManager::compile)
      .def("stats",
           [](P::StageManager &mgr) {
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
      .def(
          "run",
          [](P::StageManager &mgr, int threads, std::optional<double> flush_timeout) {
            std::optional<std::chrono::milliseconds> previous_timeout;
            if (flush_timeout.has_value()) {
              const auto override_timeout = seconds_to_milliseconds(*flush_timeout);
              previous_timeout            = mgr.get_flush_timeout();
              mgr.set_flush_timeout(override_timeout);
            }
            struct FlushTimeoutGuard {
              P::StageManager &mgr;
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
          py::arg("threads")       = 4,
          py::arg("flush_timeout") = py::none())
      .def("last_run_report",
           [](P::StageManager &mgr) -> py::object {
             const auto &opt = mgr.last_report();
             if (!opt.has_value()) return py::none();
             const auto &report = *opt;
             py::dict d;
             d["total_seconds"]             = report.total_seconds;
             d["cpu_seconds"]               = report.cpu_seconds;
             d["io_seconds"]                = report.io_seconds;
             d["ingest_seconds"]            = report.ingest_seconds;
             d["prepare_seconds"]           = report.prepare_seconds;
             d["flush_seconds"]             = report.flush_seconds;
             d["setup_seconds"]             = report.setup_seconds;
             d["compute_seconds"]           = report.compute_seconds;
             d["items_total"]               = report.items_total;
             d["items_processed"]           = report.items_processed;
             d["items_skipped"]             = report.items_skipped;
             d["stage_count"]               = report.stage_count;
             d["threads_requested"]         = report.threads_requested;
             d["threads_used"]              = report.threads_used;
             d["all_thread_safe"]           = report.all_thread_safe;
             d["run_token"]                 = report.run_token;
             d["metrics_enabled"]           = report.metrics_enabled;
             d["peak_inflight_items"]       = report.peak_inflight_items;
             d["average_queue_depth"]       = report.average_queue_depth;
             d["permit_wait_events"]        = report.permit_wait_events;
             d["permit_wait_total_seconds"] = report.permit_wait_total_seconds;
             d["permit_wait_min_seconds"]   = report.permit_wait_min_seconds;
             d["permit_wait_max_seconds"]   = report.permit_wait_max_seconds;
             d["permit_wait_avg_seconds"]   = report.permit_wait_avg_seconds;
             py::list stage_breakdown;
             for (const auto &st : report.stage_breakdown) {
               py::dict stage_entry;
               stage_entry["label"]           = st.label;
               stage_entry["setup_seconds"]   = st.setup_seconds;
               stage_entry["compute_seconds"] = st.compute_seconds;
               stage_breakdown.append(std::move(stage_entry));
             }
             d["stage_breakdown"]          = std::move(stage_breakdown);
             d["mux_sink_count"]           = report.mux_sink_count;
             d["mux_enqueued_msgs"]        = report.mux_enqueued_msgs;
             d["mux_enqueued_bytes"]       = report.mux_enqueued_bytes;
             d["mux_written_msgs"]         = report.mux_written_msgs;
             d["mux_written_bytes"]        = report.mux_written_bytes;
             d["mux_stall_ns"]             = report.mux_stall_ns;
             d["mux_drops"]                = report.mux_drops;
             d["mux_queue_depth_peak"]     = report.mux_queue_depth_peak;
             d["mux_queue_bytes_peak"]     = report.mux_queue_bytes_peak;
             d["mux_active_writers_total"] = report.mux_active_writers_total;
             d["mux_active_writers_peak"]  = report.mux_active_writers_peak;
             return d;
           })

      .def("set_auto_builtins", [](P::StageManager &mgr, bool on) { mgr.set_auto_builtins(on); })
      .def("get_auto_builtins", [](P::StageManager &mgr) { return mgr.get_auto_builtins(); })
      .def(
          "set_reporting_level",
          [](P::StageManager &mgr, P::StageManager::ReportingLevel level) { mgr.set_reporting_level(level); })
      .def("get_reporting_level", [](P::StageManager &mgr) { return mgr.get_reporting_level(); })

      .def(
          "set_flush_timeout",
          [](P::StageManager &mgr, double timeout_seconds) {
            mgr.set_flush_timeout(seconds_to_milliseconds(timeout_seconds));
          },
          py::arg("timeout"))
      .def("get_flush_timeout",
           [](P::StageManager &mgr) { return milliseconds_to_seconds(mgr.get_flush_timeout()); })

      // Parameter access for builtins
      .def("get_system_params",
           [](P::StageManager &mgr) -> py::dict {
             const auto &params = mgr.get_system_params();
             py::dict result;
             result["is_model"] = params.is_model;
             return result;
           })
      .def("set_system_params",
           [](P::StageManager &mgr, py::dict d) {
             auto &params = mgr.get_system_params();
             if (d.contains("is_model")) params.is_model = d["is_model"].cast<bool>();
             mgr.invalidate_compilation();
           })
      .def("get_topology_params",
           [](P::StageManager &mgr) -> py::dict {
             const auto &params = mgr.get_topology_params();
             py::dict result;
             result["flags"]              = static_cast<int>(params.flags);
             result["atom_typing_method"] = static_cast<int>(params.atom_typing_method);
             return result;
           })
      .def("set_topology_params",
           [](P::StageManager &mgr, py::dict d) {
             auto &params = mgr.get_topology_params();
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
                 params.atom_typing_method = static_cast<AtomTypingMethod>(
                     d["atom_typing_method"].cast<int>());
               } else {
                 params.atom_typing_method = d["atom_typing_method"].cast<AtomTypingMethod>();
               }
             }
             mgr.invalidate_compilation();
           })

      // Graph introspection
      .def("describe_graph", &P::StageManager::describe_graph);
}

} // namespace lahuta::bindings

#endif // LAHUTA_BINDINGS_STAGE_MANAGER_HPP
