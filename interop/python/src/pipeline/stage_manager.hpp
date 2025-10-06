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
#include "pipeline/dynamic/sources.hpp"
#include "pipeline/python_task.hpp"

// clang-format off
namespace py = pybind11;
namespace lahuta::bindings {
using namespace lahuta::sources;
using namespace lahuta::pipeline::dynamic;

namespace {

class PyNMRDescriptor final : public IDescriptor {
public:
  explicit PyNMRDescriptor(std::vector<std::string> paths)
      : paths_(std::move(paths)) {}

  std::optional<IngestDescriptor> next() override {
    if (index_ >= paths_.size()) return std::nullopt;
    IngestDescriptor desc;
    desc.id = paths_[index_];
    desc.origin = NMRRef{paths_[index_]};
    ++index_;
    return desc;
  }

  void reset() override { index_ = 0; }

private:
  std::vector<std::string> paths_;
  std::size_t index_ = 0;
};

class PyTrajectoryDescriptor final : public IDescriptor {
public:
  struct Entry {
    std::string id;
    std::string structure;
    std::vector<std::string> xtcs;
  };

  explicit PyTrajectoryDescriptor(std::vector<Entry> entries)
      : entries_(std::move(entries)) {}

  std::optional<IngestDescriptor> next() override {
    if (index_ >= entries_.size()) return std::nullopt;
    const auto& entry = entries_[index_++];
    IngestDescriptor desc;
    desc.id = entry.id.empty() ? entry.structure : entry.id;
    MDRef ref;
    ref.path = entry.structure;
    ref.xtc_paths = entry.xtcs;
    desc.origin = std::move(ref);
    return desc;
  }

  void reset() override { index_ = 0; }

private:
  std::vector<Entry> entries_;
  std::size_t index_ = 0;
};

inline std::vector<PyTrajectoryDescriptor::Entry>
parse_md_trajectories(py::handle obj) {
  if (!obj || obj.is_none()) {
    return {};
  }
  if (!PySequence_Check(obj.ptr())) {
    throw std::invalid_argument("from_md_trajectories expects a sequence");
  }
  std::vector<PyTrajectoryDescriptor::Entry> entries;
  py::sequence seq = py::reinterpret_borrow<py::sequence>(obj);
  const std::size_t n = py::len(seq);
  entries.reserve(n);
  for (auto item : seq) {
    std::string structure;
    std::vector<std::string> xtcs;
    std::string id;

    if (py::isinstance<py::dict>(item)) {
      py::dict d = item.cast<py::dict>();
      if (d.contains("structure")) {
        structure = d[py::str("structure")].cast<std::string>();
      } else if (d.contains("gro")) {
        structure = d[py::str("gro")].cast<std::string>();
      } else {
        throw std::invalid_argument("from_md_trajectories expects 'structure' or 'gro' in each mapping");
      }
      py::object xtc_obj;
      if (d.contains("xtc")) xtc_obj = d[py::str("xtc")];
      else if (d.contains("xtcs")) xtc_obj = d[py::str("xtcs")];
      else throw std::invalid_argument("from_md_trajectories expects 'xtc' or 'xtcs' in each mapping");
      if (d.contains("id")) id = d[py::str("id")].cast<std::string>();
      else if (d.contains("session_id")) id = d[py::str("session_id")].cast<std::string>();
      xtcs = [&]() {
        std::vector<std::string> out;
        if (py::isinstance<py::str>(xtc_obj) || py::isinstance<py::bytes>(xtc_obj)) {
          out.push_back(xtc_obj.cast<std::string>());
        } else {
          py::sequence xtc_seq = xtc_obj.cast<py::sequence>();
          out.reserve(py::len(xtc_seq));
          for (auto v : xtc_seq) out.push_back(v.cast<std::string>());
        }
        return out;
      }();
    } else {
      py::sequence tup = py::reinterpret_borrow<py::sequence>(item);
      if (py::len(tup) < 2) {
        throw std::invalid_argument("from_md_trajectories expects (structure, xtc[, id]) tuples");
      }
      structure = tup[0].cast<std::string>();
      py::object xtc_obj = tup[1];
      if (py::len(tup) > 2) id = tup[2].cast<std::string>();
      if (py::isinstance<py::str>(xtc_obj) || py::isinstance<py::bytes>(xtc_obj)) {
        xtcs.push_back(xtc_obj.cast<std::string>());
      } else {
        py::sequence xtc_seq = xtc_obj.cast<py::sequence>();
        xtcs.reserve(py::len(xtc_seq));
        for (auto v : xtc_seq) xtcs.push_back(v.cast<std::string>());
      }
    }

    if (structure.empty()) {
      throw std::invalid_argument("from_md_trajectories requires non-empty structure paths");
    }
    if (xtcs.empty()) {
      throw std::invalid_argument("from_md_trajectories requires at least one XTC path per entry");
    }

    entries.push_back(PyTrajectoryDescriptor::Entry{
        id.empty() ? structure : id,
        std::move(structure),
        std::move(xtcs)
    });
  }
  return entries;
}

} // namespace

inline void bind_stage_manager(py::module_ &md) {
  py::class_<ITask, std::shared_ptr<ITask>> task(md, "Task");

  // Bind reusable C++ task that serializes models to ModelRecord payloads
  py::class_<analysis::system::ModelPackTask, ITask, std::shared_ptr<analysis::system::ModelPackTask>>(md, "ModelPackTask")
    .def(py::init<std::string>(), py::arg("channel") = std::string("db"));

  py::class_<StageManager, std::shared_ptr<StageManager>>(md, "StageManager")
    // sources
    .def_static("from_directory", [](const std::string &path, const std::string &ext, bool recursive, std::size_t batch) {
          auto src = lahuta::pipeline::dynamic::sources_factory::from_directory(path, ext, recursive, batch);
          return std::make_shared<StageManager>(std::move(src));
        },
        py::arg("path"), py::arg("ext") = std::string(""), py::arg("recursive") = true, py::arg("batch") = 200)
    .def_static("from_nmr_files", [](std::vector<std::string> files) {
                  auto src = std::make_unique<PyNMRDescriptor>(std::move(files));
                  return std::make_shared<StageManager>(std::move(src));
                },
                py::arg("paths"),
                py::doc(R"doc(Create a StageManager that streams multi-model NMR structures.)doc"))
    .def_static("from_files", [](std::vector<std::string> files) {
                  auto src = lahuta::pipeline::dynamic::sources_factory::from_vector(std::move(files));
                  return std::make_shared<StageManager>(std::move(src));
                })
    .def_static("from_filelist", [](const std::string &list_path) {
                  auto src = lahuta::pipeline::dynamic::sources_factory::from_filelist(list_path);
                  return std::make_shared<StageManager>(std::move(src));
                })
    .def_static("from_md_trajectories", [](py::object trajectories) {
                  auto entries = parse_md_trajectories(trajectories);
                  auto src = std::make_unique<PyTrajectoryDescriptor>(std::move(entries));
                  return std::make_shared<StageManager>(std::move(src));
                },
                py::arg("trajectories"),
                py::doc(R"doc(Create a StageManager that streams MD trajectories (structure + XTC).)doc"))
    .def_static("from_database", [](const std::string &db_path, std::size_t batch) {
                  auto src = lahuta::pipeline::dynamic::sources_factory::from_lmdb(db_path, std::string{}, batch);
                  return std::make_shared<StageManager>(std::move(src));
                },
                py::arg("path"), py::arg("batch") = 1024)
    .def_static("from_database_handle", [](std::shared_ptr<LMDBDatabase> db, std::size_t batch) {
                  auto src = lahuta::pipeline::dynamic::sources_factory::from_lmdb(std::move(db), std::string{}, batch);
                  return std::make_shared<StageManager>(std::move(src));
                },
                py::arg("db"), py::arg("batch") = 1024)
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
