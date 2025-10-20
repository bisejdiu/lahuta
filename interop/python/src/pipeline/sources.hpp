#ifndef LAHUTA_BINDINGS_PIPELINE_SOURCES_HPP
#define LAHUTA_BINDINGS_PIPELINE_SOURCES_HPP

#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "db/db.hpp"
#include "sources/adapters/directory.hpp"
#include "sources/adapters/file_list.hpp"
#include "sources/adapters/lmdb.hpp"
#include "sources/adapters/vector.hpp"
#include "sources/descriptor.hpp"
#include "sources/directory.hpp"

namespace py = pybind11;

// clang-format off
namespace lahuta::bindings {

class PyNMRDescriptor final : public sources::IDescriptor {
public:
  explicit PyNMRDescriptor(std::vector<std::string> paths) : paths_(std::move(paths)) {}

  std::optional<IngestDescriptor> next() override {
    if (index_ >= paths_.size()) return std::nullopt;
    IngestDescriptor desc;
    desc.id     = paths_[index_];
    desc.origin = NMRRef{paths_[index_]};
    ++index_;
    return desc;
  }

  void reset() override { index_ = 0; }

private:
  std::vector<std::string> paths_;
  std::size_t index_ = 0;
};

class PyTrajectoryDescriptor final : public sources::IDescriptor {
public:
  struct Entry {
    std::string id;
    std::string structure;
    std::vector<std::string> xtcs;
  };

  explicit PyTrajectoryDescriptor(std::vector<Entry> entries) : entries_(std::move(entries)) {}

  std::optional<IngestDescriptor> next() override {
    if (index_ >= entries_.size()) return std::nullopt;
    const auto &entry = entries_[index_++];
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

inline std::vector<PyTrajectoryDescriptor::Entry> parse_md_trajectories(py::handle obj) {

  if (!obj || obj.is_none()) return {};
  if (!PySequence_Check(obj.ptr())) throw std::invalid_argument("MdTrajectoriesSource expects a sequence");

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
        throw std::invalid_argument("MdTrajectoriesSource expects 'structure' or 'gro' in each mapping");
      }

      py::object xtc_obj;
      if      (d.contains("xtc"))  xtc_obj = d[py::str("xtc")];
      else if (d.contains("xtcs")) xtc_obj = d[py::str("xtcs")];
      else throw std::invalid_argument("MdTrajectoriesSource expects 'xtc' or 'xtcs' in each mapping");

      if (d.contains("id")) {
        id = d[py::str("id")].cast<std::string>();
      }
      else if (d.contains("session_id")) {
        id = d[py::str("session_id")].cast<std::string>();
      }
      if (py::isinstance<py::str>(xtc_obj) || py::isinstance<py::bytes>(xtc_obj)) {
        xtcs.push_back(xtc_obj.cast<std::string>());
      } else {
        py::sequence xtc_seq = xtc_obj.cast<py::sequence>();
        xtcs.reserve(py::len(xtc_seq));
        for (auto v : xtc_seq) {
          xtcs.push_back(v.cast<std::string>());
        }
      }
    } else {
      py::sequence tup = py::reinterpret_borrow<py::sequence>(item);
      if (py::len(tup) < 2) throw std::invalid_argument("MdTrajectoriesSource expects (structure, xtc[, id]) tuples");

      structure = tup[0].cast<std::string>();
      py::object xtc_obj = tup[1];
      if (py::len(tup) > 2) id = tup[2].cast<std::string>();

      if (py::isinstance<py::str>(xtc_obj) || py::isinstance<py::bytes>(xtc_obj)) {
        xtcs.push_back(xtc_obj.cast<std::string>());
      } else {
        py::sequence xtc_seq = xtc_obj.cast<py::sequence>();
        xtcs.reserve(py::len(xtc_seq));
        for (auto v : xtc_seq)
          xtcs.push_back(v.cast<std::string>());
      }
    }

    if (structure.empty()) throw std::invalid_argument("MdTrajectoriesSource requires non-empty structure paths");
    if (xtcs.empty())      throw std::invalid_argument("MdTrajectoriesSource requires at least one XTC path per entry");

    entries.push_back(PyTrajectoryDescriptor::Entry{id.empty() ? structure : id, std::move(structure), std::move(xtcs)});
  }
  return entries;
}

class HandleLmdbDescriptor final : public sources::IDescriptor {
public:
  HandleLmdbDescriptor(std::shared_ptr<LMDBDatabase> db, const std::string &db_name, std::size_t batch)
      : inner_(std::move(db), db_name, batch) {}

  std::optional<IngestDescriptor> next() override { return inner_.next(); }
  void reset() override { inner_.reset(); }

private:
  sources::LMDBAdapter inner_;
};

inline void bind_sources(py::module_ &md) {
  py::module_ ms = md.def_submodule("sources", "Pipeline sources");

  // Base type w unique ownership so it can be moved into StageManager
  py::class_<sources::IDescriptor, std::shared_ptr<sources::IDescriptor>> _(ms, "Source");

  py::class_<sources::DirectoryAdapter, sources::IDescriptor, std::shared_ptr<sources::DirectoryAdapter>>(ms, "DirectorySource")
      .def(py::init<const std::string &, const std::string &, bool, std::size_t>(),
           py::arg("path"), py::arg("ext") = std::string(""), py::arg("recursive") = true, py::arg("batch") = 200)
      .def(py::init([](const std::string &path, std::vector<std::string> extensions, bool recursive, std::size_t batch) {
             return std::make_shared<sources::DirectoryAdapter>(path, std::move(extensions), recursive, batch);
           }),
           py::arg("path"), py::arg("extensions"), py::arg("recursive") = true, py::arg("batch") = 200);

  py::class_<sources::VectorAdapter, sources::IDescriptor, std::shared_ptr<sources::VectorAdapter>>(ms, "FileSource")
      .def(py::init<std::vector<std::string>>(), py::arg("files"));

  py::class_<sources::FileListAdapter, sources::IDescriptor, std::shared_ptr<sources::FileListAdapter>>(ms, "FileListSource")
      .def(py::init<const std::string &>(), py::arg("path"));

  // LMDB sources with env-path variant
  py::class_<sources::LMDBAdapter, sources::IDescriptor, std::shared_ptr<sources::LMDBAdapter>>(ms, "DatabaseSource")
      .def(py::init<const std::string &, const std::string &, std::size_t>(),
           py::arg("path"), py::arg("database") = std::string(""), py::arg("batch") = 1024);

  // LMDB sources with handle variant
  py::class_<HandleLmdbDescriptor, sources::IDescriptor, std::shared_ptr<HandleLmdbDescriptor>>(ms, "DatabaseHandleSource")
      .def(py::init<std::shared_ptr<LMDBDatabase>, const std::string &, std::size_t>(),
           py::arg("db"), py::arg("database") = std::string(""), py::arg("batch") = 1024);

  // NMR and MD streaming sources
  py::class_<PyNMRDescriptor, sources::IDescriptor, std::shared_ptr<PyNMRDescriptor>>(ms, "NmrSource")
      .def(py::init<std::vector<std::string>>(), py::arg("files"),
          py::doc(R"doc(Create a source that streams multi-model NMR structures.)doc"));

  py::class_<PyTrajectoryDescriptor, sources::IDescriptor, std::shared_ptr<PyTrajectoryDescriptor>>(ms, "MdTrajectoriesSource")
      .def(py::init([](py::object trajectories) {
            return std::make_shared<PyTrajectoryDescriptor>(parse_md_trajectories(trajectories));
          }),
          py::arg("trajectories"),
          py::doc(R"doc(Create a source that streams MD trajectories (structure + XTC).)doc"));
}

} // namespace lahuta::bindings

#endif // LAHUTA_BINDINGS_PIPELINE_SOURCES_HPP
