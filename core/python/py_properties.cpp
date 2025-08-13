#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/RDKitBase.h>

#include "array.hpp"
#include "contacts.hpp"
#include "lahuta.hpp"
#include "luni_processor.hpp"
#include "properties/analyzer.hpp"
#include "properties/query.hpp"
#include "properties/types.hpp"
#include "py_properties.hpp"

// clang-format off
using namespace lahuta;
namespace py = pybind11;

using LuniResultType    = PropertyResult<std::vector<int>, std::vector<std::string>, std::vector<RDGeom::Point3D>>;
using LuniAnalyzer      = PropertyAnalyzer<Luni>;
using LuniFileProcessor = FileProcessor<Luni, LuniAnalyzer>;

using LuniIdentityAnalyzer  = IdentityAnalyzer<Luni>; // FIX: a temporary solution
using LuniFileProcessor2    = FileProcessor<Luni, LuniIdentityAnalyzer>;

using GetterFunc = py::array (*)(const LuniResultType &);

struct GetterEntry {
  PropertyKey key;
  GetterFunc func;
};

py::array get_names (const LuniResultType &result) {
  return string_array(result.get<PropertyKey::Names>());
}

py::array get_indices  (const LuniResultType &result) {
  return int_array(result.get<PropertyKey::Indices>());
}

py::array get_elements (const LuniResultType &result) {
  return string_array(result.get<PropertyKey::Elements>());
}

py::array get_positions(const LuniResultType &result) {
  return coordinates(result.get<PropertyKey::Positions>());
}

constexpr std::array<GetterEntry, 4> GetterEntries = {{
    {PropertyKey::Names,      get_names},
    {PropertyKey::Indices,    get_indices},
    {PropertyKey::Elements,   get_elements},
    {PropertyKey::Positions,  get_positions},
}};

auto get_property = [](const LuniResultType &self, PropertyKey key) -> py::array {
  auto it = std::find_if(GetterEntries.begin(), GetterEntries.end(), [key](const GetterEntry& entry) { return entry.key == key; });
  if (it == GetterEntries.end()) {
    throw py::key_error("Unknown property key");
  }
  return it->func(self);
};

void bind_properties(py::module &m) {

  py::enum_<PropertyKey> PropertyKey_(m, "PropertyKey");

  py::class_<PropertyQuery<Luni>>     PropertyQuery_    (m, "PropertyQueryLuni");
  py::class_<PropertyAnalyzer<Luni>>  PropertyAnalyzer_ (m, "PropertyAnalyzerLuni");
  py::class_<IdentityAnalyzer<Luni>>  IdentityAnalyzer_ (m, "IdentityAnalyzerLuni");
  py::class_<LuniResultType>          LuniResultType_   (m, "LuniPropertyResult");
  py::class_<LuniFileProcessor>       LuniFileProcessor_(m, "LuniFileProcessor");


  PropertyKey_
      .value("Names",     PropertyKey::Names)
      .value("Indices",   PropertyKey::Indices)
      .value("Elements",  PropertyKey::Elements)
      .value("Positions", PropertyKey::Positions);

  PropertyQuery_
      .def(py::init<>())
      .def("select", [](PropertyQuery<Luni> &self, PropertyKey key) -> PropertyQuery<Luni> & { return self.select(key); })
      .def("select", [](PropertyQuery<Luni> &self, const std::vector<PropertyKey> &keys) -> PropertyQuery<Luni> & {
            for (auto &key : keys) {
              self.select(key);
            }
            return self;
          })
      .def("properties", &PropertyQuery<Luni>::properties, py::return_value_policy::reference_internal);

  PropertyAnalyzer_
      .def(py::init<PropertyQuery<Luni>>());

  IdentityAnalyzer_
      .def(py::init<>());

  LuniResultType_
    .def("get", get_property)
    .def("has_property", [](const LuniResultType &self, PropertyKey key) {
        auto it = std::find_if(GetterEntries.begin(), GetterEntries.end(), [key](const GetterEntry& entry) { return entry.key == key; });
        return it != GetterEntries.end();
    })
    .def("__getitem__", get_property)
    .def("__contains__", [](const LuniResultType &self, PropertyKey key) -> bool {
        auto it = std::find_if(GetterEntries.begin(), GetterEntries.end(),
                               [key](const GetterEntry& entry) { return entry.key == key; });
        return it != GetterEntries.end();
    });

  LuniFileProcessor_
      .def(py::init<int, LuniAnalyzer, bool>(), py::arg("concurrency"), py::arg("analyzer"), py::arg("use_progress_bar") = true)
      .def("process_files",       &LuniFileProcessor::process_files)
      .def("wait_for_completion", &LuniFileProcessor::wait_for_completion)

      .def("get_result", [](LuniFileProcessor &self, const std::string &file_name) -> py::object {
            auto *res_ptr = self.get_result(file_name);
            if (!res_ptr) {
              return py::none();
            }

            return py::cast(res_ptr, py::return_value_policy::reference);
          })

      .def("get_all_results", [](LuniFileProcessor &self) {
        auto raw_map = self.get_all_results(); // returns std::unordered_map<std::string, ResultType*>
        py::dict py_map;
        for (auto &kv : raw_map) {
          const std::string &fname = kv.first;
          auto *res_ptr = kv.second;
          if (res_ptr) {
            py_map[py::str(fname.c_str())] = py::cast(res_ptr, py::return_value_policy::reference);
          } else {
            py_map[py::str(fname.c_str())] = py::none();
          }
        }
        return py_map;
      });

  m.def("process_files", [](const std::vector<std::string> &files, const std::vector<PropertyKey> &property_keys, int concurrency, bool use_spinner) -> py::dict {
        PropertyQuery<Luni> query;
        for (const auto &key : property_keys) {
          query.select(key);
        }

        PropertyAnalyzer<Luni> analyzer(query);
        LuniFileProcessor processor(concurrency, analyzer, use_spinner);

        processor.process_files(files);
        processor.wait_for_completion();

        py::dict results;
        for (auto &f : files) {
          auto *res_ptr = processor.get_result(f);
          if (!res_ptr) {
            results[py::str(f)] = py::none();
            continue;
          }
          // Move the result out of the processor-owned unique_ptr
          auto moved_result = std::move(*res_ptr);
          results[py::str(f)] = py::cast(std::move(moved_result), py::return_value_policy::move);
        }
        return results;
      }, py::arg("files"), py::arg("property_keys"), py::arg("n_jobs"), py::arg("use_spinner") = true);

  m.def("process_files", [](const std::vector<std::string> &files, int concurrency, bool use_spinner) -> py::dict {

        IdentityAnalyzer<Luni> analyzer;
        LuniFileProcessor2 processor(concurrency, analyzer, use_spinner);

        processor.process_files(files);
        processor.wait_for_completion();

        py::dict results;
        for (auto &f : files) {
          auto *res_ptr = processor.get_result(f);
          if (!res_ptr) {
            results[py::str(f)] = py::none();
            continue;
          }
          // Move the result out of the processor-owned unique_ptr
          auto moved_result = std::move(*res_ptr);
          results[py::str(f)] = py::cast(std::move(moved_result), py::return_value_policy::move);
        }
        return results;
      }, py::arg("files"), py::arg("n_jobs"), py::arg("use_spinner") = true);
}
