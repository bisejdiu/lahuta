#include <memory>
#include <stdexcept>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "analysis/system/records.hpp"
#include "bindings.hpp"
#include "db/db.hpp"
#include "db/reader.hpp"
#include "lahuta.hpp"
#include "models/topology.hpp"

namespace py = pybind11;
namespace lahuta::bindings {

// clang-format off
void bind_db(py::module_ &m) {
  py::module_ mdb = m.def_submodule("db", "LMDB integration");

  py::class_<LMDBDatabase, std::shared_ptr<LMDBDatabase>>(mdb, "Database")
    .def(py::init<const std::string&, std::size_t>(), py::arg("path"), py::arg("max_size_gb") = 500)
    .def("keys", [](LMDBDatabase &db) {
          std::vector<std::string> out;
          db.for_each_key([&](const std::string &k) { out.push_back(k); });
          return out;
        },
        "Return all keys in the database (eager list)")
    .def("get_raw", [](LMDBDatabase &db, const std::string &key) {
          LMDBReader r(db.get_env(), db.get_dbi());
          std::string_view raw;
          if (!r.fetch(key, raw)) throw std::runtime_error("Key not found: " + key);
          return py::bytes(raw.data(), raw.size());
        },
        py::arg("key"))
    .def("get_model", [](LMDBDatabase &db, const std::string &key) {
          LMDBReader r(db.get_env(), db.get_dbi());
          std::string_view raw;
          if (!r.fetch(key, raw)) throw std::runtime_error("Key not found: " + key);
          using serialization::Serializer;
          using lahuta::analysis::system::ModelRecord;
          auto rec = Serializer<fmt::binary, ModelRecord>::deserialize(raw.data(), raw.size());
          auto mol = std::make_shared<RDKit::RWMol>();
          if (!lahuta::build_model_topology(mol, rec.data, ModelTopologyMethod::CSR)) {
            throw std::runtime_error("Failed to build model topology");
          }
          auto sys = Luni::create(mol);
          return std::make_shared<Luni>(std::move(sys));
        },
        py::arg("key"),
        "Create and return a LahutaSystem built from the stored model");
}

} // namespace lahuta::bindings
