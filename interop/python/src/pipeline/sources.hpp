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
#include "pipeline/dynamic/sources.hpp"
#include "sources/adapters/directory.hpp"
#include "sources/descriptor.hpp"

namespace py = pybind11;

namespace lahuta::bindings {

using namespace lahuta::pipeline::dynamic;

class PySource {
public:
  explicit PySource(std::unique_ptr<sources::IDescriptor> descriptor)
      : descriptor_(std::move(descriptor)) {
    if (!descriptor_) {
      throw std::invalid_argument("Pipeline source requires a descriptor");
    }
  }
  virtual ~PySource() = default;

  std::unique_ptr<sources::IDescriptor> release() {
    if (!descriptor_) {
      throw std::runtime_error("Pipeline source has already been consumed");
    }
    return std::move(descriptor_);
  }

  bool available() const noexcept { return static_cast<bool>(descriptor_); }

protected:
  void reset_descriptor(std::unique_ptr<sources::IDescriptor> descriptor) {
    descriptor_ = std::move(descriptor);
  }

private:
  std::unique_ptr<sources::IDescriptor> descriptor_;
};

class DirectorySource final : public PySource {
public:
  DirectorySource(const std::string &path, const std::string &ext, bool recursive, std::size_t batch)
      : PySource(sources_factory::from_directory(path, ext, recursive, batch)) {}
};

class FilesSource final : public PySource {
public:
  explicit FilesSource(std::vector<std::string> files)
      : PySource(sources_factory::from_vector(std::move(files))) {}
};

class FileListSource final : public PySource {
public:
  explicit FileListSource(const std::string &list_path)
      : PySource(sources_factory::from_filelist(list_path)) {}
};

class DatabaseSource final : public PySource {
public:
  DatabaseSource(const std::string &env_path, const std::string &db_name, std::size_t batch)
      : PySource(sources_factory::from_lmdb(env_path, db_name, batch)) {}
};

class DatabaseHandleSource final : public PySource {
public:
  DatabaseHandleSource(std::shared_ptr<LMDBDatabase> db, const std::string &db_name, std::size_t batch)
      : PySource(sources_factory::from_lmdb(std::move(db), db_name, batch)) {}
};

class PyNMRDescriptor final : public sources::IDescriptor {
public:
  explicit PyNMRDescriptor(std::vector<std::string> paths)
      : paths_(std::move(paths)) {}

  std::optional<lahuta::IngestDescriptor> next() override {
    if (index_ >= paths_.size()) return std::nullopt;
    lahuta::IngestDescriptor desc;
    desc.id = paths_[index_];
    desc.origin = lahuta::NMRRef{paths_[index_]};
    ++index_;
    return desc;
  }

  void reset() override { index_ = 0; }

private:
  std::vector<std::string> paths_;
  std::size_t index_ = 0;
};

class NmrSource final : public PySource {
public:
  explicit NmrSource(std::vector<std::string> files)
      : PySource(std::make_unique<PyNMRDescriptor>(std::move(files))) {}
};

class PyTrajectoryDescriptor final : public sources::IDescriptor {
public:
  struct Entry {
    std::string id;
    std::string structure;
    std::vector<std::string> xtcs;
  };

  explicit PyTrajectoryDescriptor(std::vector<Entry> entries)
      : entries_(std::move(entries)) {}

  std::optional<lahuta::IngestDescriptor> next() override {
    if (index_ >= entries_.size()) return std::nullopt;
    const auto &entry = entries_[index_++];
    lahuta::IngestDescriptor desc;
    desc.id = entry.id.empty() ? entry.structure : entry.id;
    lahuta::MDRef ref;
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
  if (!obj || obj.is_none()) {
    return {};
  }
  if (!PySequence_Check(obj.ptr())) {
    throw std::invalid_argument("MdTrajectoriesSource expects a sequence");
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
        throw std::invalid_argument("MdTrajectoriesSource expects 'structure' or 'gro' in each mapping");
      }
      py::object xtc_obj;
      if (d.contains("xtc")) xtc_obj = d[py::str("xtc")];
      else if (d.contains("xtcs")) xtc_obj = d[py::str("xtcs")];
      else throw std::invalid_argument("MdTrajectoriesSource expects 'xtc' or 'xtcs' in each mapping");
      if (d.contains("id")) id = d[py::str("id")].cast<std::string>();
      else if (d.contains("session_id")) id = d[py::str("session_id")].cast<std::string>();
      if (py::isinstance<py::str>(xtc_obj) || py::isinstance<py::bytes>(xtc_obj)) {
        xtcs.push_back(xtc_obj.cast<std::string>());
      } else {
        py::sequence xtc_seq = xtc_obj.cast<py::sequence>();
        xtcs.reserve(py::len(xtc_seq));
        for (auto v : xtc_seq) xtcs.push_back(v.cast<std::string>());
      }
    } else {
      py::sequence tup = py::reinterpret_borrow<py::sequence>(item);
      if (py::len(tup) < 2) {
        throw std::invalid_argument("MdTrajectoriesSource expects (structure, xtc[, id]) tuples");
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
      throw std::invalid_argument("MdTrajectoriesSource requires non-empty structure paths");
    }
    if (xtcs.empty()) {
      throw std::invalid_argument("MdTrajectoriesSource requires at least one XTC path per entry");
    }

    entries.push_back(PyTrajectoryDescriptor::Entry{
        id.empty() ? structure : id,
        std::move(structure),
        std::move(xtcs)
    });
  }
  return entries;
}

class MdTrajectoriesSource final : public PySource {
public:
  explicit MdTrajectoriesSource(py::object trajectories)
      : PySource(std::make_unique<PyTrajectoryDescriptor>(parse_md_trajectories(trajectories))) {}
};

inline void bind_sources(py::module_ &md) {
  py::module_ ms = md.def_submodule("sources", "Pipeline sources");

  py::class_<PySource, std::shared_ptr<PySource>>(ms, "Source")
      .def_property_readonly("available", &PySource::available,
                             "Whether the source descriptor is still available for consumption.");

  py::class_<DirectorySource, PySource, std::shared_ptr<DirectorySource>>(ms, "DirectorySource")
      .def(py::init<const std::string &, const std::string &, bool, std::size_t>(),
           py::arg("path"), py::arg("ext") = std::string(""), py::arg("recursive") = true, py::arg("batch") = 200);

  py::class_<FilesSource, PySource, std::shared_ptr<FilesSource>>(ms, "FilesSource")
      .def(py::init<std::vector<std::string>>(), py::arg("files"));

  py::class_<FileListSource, PySource, std::shared_ptr<FileListSource>>(ms, "FileListSource")
      .def(py::init<const std::string &>(), py::arg("path"));

  py::class_<DatabaseSource, PySource, std::shared_ptr<DatabaseSource>>(ms, "DatabaseSource")
      .def(py::init<const std::string &, const std::string &, std::size_t>(),
           py::arg("path"), py::arg("database") = std::string(""), py::arg("batch") = 1024);

  py::class_<DatabaseHandleSource, PySource, std::shared_ptr<DatabaseHandleSource>>(ms, "DatabaseHandleSource")
      .def(py::init<std::shared_ptr<LMDBDatabase>, const std::string &, std::size_t>(),
           py::arg("db"), py::arg("database") = std::string(""), py::arg("batch") = 1024);

  py::class_<NmrSource, PySource, std::shared_ptr<NmrSource>>(ms, "NmrSource")
      .def(py::init<std::vector<std::string>>(), py::arg("files"),
           py::doc(R"doc(Create a source that streams multi-model NMR structures.)doc"));

  py::class_<MdTrajectoriesSource, PySource, std::shared_ptr<MdTrajectoriesSource>>(ms, "MdTrajectoriesSource")
      .def(py::init<py::object>(), py::arg("trajectories"),
           py::doc(R"doc(Create a source that streams MD trajectories (structure + XTC).)doc"));
}

} // namespace lahuta::bindings

#endif // LAHUTA_BINDINGS_PIPELINE_SOURCES_HPP
