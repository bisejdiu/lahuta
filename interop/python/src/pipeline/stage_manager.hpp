#ifndef LAHUTA_BINDINGS_STAGE_MANAGER_HPP
#define LAHUTA_BINDINGS_STAGE_MANAGER_HPP

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "analysis/contacts/computation.hpp"
#include "analysis/contacts/provider.hpp"
#include "pipeline/compute/parameters.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/python_task.hpp"

// clang-format off
namespace py = pybind11;
namespace lahuta::bindings {
using namespace lahuta::sources;

inline void bind_stage_manager(py::module_ &md) {
  py::class_<ITask, std::shared_ptr<ITask>> task(md, "Task");

  py::class_<StageManager, std::shared_ptr<StageManager>>(md, "StageManager")
    // sources
    .def_static("from_directory", [](const std::string &path, const std::string &ext, bool recursive, std::size_t batch) {
          return std::make_shared<StageManager>(StageManager::SourceVariant(
              std::in_place_type<DirectorySource>, path, ext, recursive, batch));
        },
        py::arg("path"), py::arg("ext") = std::string(""), py::arg("recursive") = true, py::arg("batch") = 200)
    .def_static("from_files", [](std::vector<std::string> files) {
                  return std::make_shared<StageManager>(
                      StageManager::SourceVariant(
                          std::in_place_type<VectorSource>,
                          std::move(files)));
                })
    .def_static("from_filelist", [](const std::string &list_path) {
                  return std::make_shared<StageManager>(
                      StageManager::SourceVariant(
                          std::in_place_type<lahuta::sources::FileListSource>,
                          list_path));
                })
    .def_static("from_database", [](const std::string &db_path, std::size_t batch) {
                  return std::make_shared<StageManager>(
                      StageManager::SourceVariant(
                          std::in_place_type<DBKeySourceHolder>,
                          db_path, batch));
                },
                py::arg("path"), py::arg("batch") = 1024)
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
                            lahuta::analysis::contacts::ContactProvider provider,
                            InteractionType interaction_type,
                            std::optional<std::string> channel,
                            const std::string &out_fmt,
                            bool thread_safe) {
          std::string ch = channel.has_value() ? *channel : name;
          bool json = true;
          if (out_fmt == "json") json = true;
          else if (out_fmt == "text") json = false;
          else throw std::invalid_argument("out_fmt must be 'json' or 'text'");

          pipeline::compute::ContactsParams p{};
          p.provider = provider;
          p.type     = interaction_type;
          p.channel  = ch;
          p.json     = json;

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
          mgr.add_task(name, deps, std::move(t), /*thread_safe_for_stage=*/true);
        },
        py::arg("name"),
        py::arg("depends") = std::vector<std::string>{},
        py::arg("fn"),
        py::arg("channel") = std::optional<std::string>{},
        py::arg("serialize") = true,
        py::arg("store") = true)

    // sinks
    .def("connect_sink", [](StageManager &mgr, const std::string &channel, std::shared_ptr<IDynamicSink> sink) {
          mgr.connect_sink(channel, std::move(sink));
        },
        py::arg("channel"), py::arg("sink"))
    .def("sorted_tasks", &StageManager::sorted_tasks)
    .def("compile",      &StageManager::compile)
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
          mgr.invalidate_compilation();
        })

    // Graph introspection
    .def("describe_graph", &StageManager::describe_graph);
}

} // namespace lahuta::bindings

#endif // LAHUTA_BINDINGS_STAGE_MANAGER_HPP
