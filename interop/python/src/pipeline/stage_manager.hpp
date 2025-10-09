#ifndef LAHUTA_BINDINGS_STAGE_MANAGER_HPP
#define LAHUTA_BINDINGS_STAGE_MANAGER_HPP

#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "analysis/contacts/computation.hpp"
#include "analysis/contacts/provider.hpp"
#include "analysis/system/model_pack_task.hpp"
#include "pipeline/compute/parameters.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/python_task.hpp"

// clang-format off
namespace py = pybind11;
namespace lahuta::bindings {
using namespace lahuta::sources;
using namespace lahuta::pipeline::dynamic;

inline void bind_stage_manager(py::module_ &md) {
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
                            const std::vector<std::string> &deps,
                            analysis::contacts::ContactProvider provider,
                            InteractionType interaction_type,
                            std::optional<std::string> channel,
                            const std::string &out_fmt,
                            bool thread_safe) {
          std::string ch = channel.has_value() ? *channel : name;

          pipeline::compute::ContactsOutputFormat format;
          if      (out_fmt == "json")   format = pipeline::compute::ContactsOutputFormat::Json;
          else if (out_fmt == "text")   format = pipeline::compute::ContactsOutputFormat::Text;
          else if (out_fmt == "binary") format = pipeline::compute::ContactsOutputFormat::Binary;
          else throw std::invalid_argument("out_fmt must be 'json', 'text', or 'binary'");

          pipeline::compute::ContactsParams p{};
          p.provider = provider;
          p.type     = interaction_type;
          p.channel  = ch;
          p.format   = format;

          mgr.add_computation(name, deps, [label = name, p]() {
            return std::make_unique<analysis::contacts::ContactsComputation>(label, p);
          }, thread_safe);
        },
        py::arg("name"), py::arg("depends") = std::vector<std::string>{},
        py::arg("provider"), py::arg("interaction_type"),
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
            d["sink_name"]     = s.sink_name;
            d["enqueued_msgs"] = s.enqueued_msgs;
            d["enqueued_bytes"]= s.enqueued_bytes;
            d["written_msgs"]  = s.written_msgs;
            d["written_bytes"] = s.written_bytes;
            d["stalled_ns"]    = s.stalled_ns;
            d["drops"]         = s.drops;
            d["queue_msgs"]    = s.queue_msgs;
            d["queue_bytes"]   = s.queue_bytes;
            out.append(std::move(d));
          }
          return out;
        })
    .def("run", [](StageManager &mgr, int threads) {
          {
            py::gil_scoped_release release;
            mgr.run(threads);
          }
          py::gil_scoped_acquire acquire;
        },
        py::arg("threads") = 8)

    // Policy control for builtins. Python turns this ON by default
    .def("set_auto_builtins", [](StageManager& mgr, bool on) { mgr.set_auto_builtins(on); })
    .def("get_auto_builtins", [](StageManager& mgr) { return mgr.get_auto_builtins(); })

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
