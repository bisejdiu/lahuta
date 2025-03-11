#include <optional>
#include <pybind11/complex.h> // FIX: is this needed?
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/RDKitBase.h>

#include "at.hpp"
#include "common.hpp"
#include "contacts.hpp"
#include "lahuta.hpp"
#include "entities.hpp"
#include "nsgrid.hpp"
#include "residues.hpp"
#include "topology.hpp"
#include "distances.hpp"
#include "distopia.h"
#include "ent.hpp"

namespace py = pybind11;
using namespace lahuta;

// clang-format off

template <typename T>
static py::array_t<T> distance(const std::vector<std::vector<T>> &points1, const std::vector<std::vector<T>> &points2) {
  const int na = static_cast<int>(points1.size());
  const int nb = static_cast<int>(points2.size());

  if (na == 0 || nb == 0) return py::array_t<T>();

  std::vector<T> a(3 * na);
  for (int i = 0; i < na; ++i) {
    a[3 * i + 0] = points1[i][0];
    a[3 * i + 1] = points1[i][1];
    a[3 * i + 2] = points1[i][2];
  }

  std::vector<T> b(3 * nb);
  for (int i = 0; i < nb; ++i) {
    b[3 * i + 0] = points2[i][0];
    b[3 * i + 1] = points2[i][1];
    b[3 * i + 2] = points2[i][2];
  }

  auto result = py::array_t<T>({na, nb});
  auto buf = result.request();
  T* ptr = static_cast<T*>(buf.ptr);

  distopia::DistanceArrayNoBox(a.data(), b.data(), na, nb, ptr);

  return result;
}

template <typename T>
static py::array_t<T> distance(const std::vector<std::vector<T>> &points) {
  return distance(points, points);
}

void bind(py::module &_lahuta) {
  py::class_<RDKit::Atom> Atom_(_lahuta, "Atom");

  py::class_<Luni>                Luni       (_lahuta, "LahutaCPP");
  py::class_<Topology>            Topology_  (_lahuta, "Topology_");
  py::class_<Residue>             Residue_   (_lahuta, "Residue");
  py::class_<Residues>            Residues_  (_lahuta, "Residues_");
  py::class_<IR>                  IntRepr_   (_lahuta, "IR");
  py::class_<NSResults>           NSResults_ (_lahuta, "NSResults_");
  py::class_<FastNS>              FastNS_    (_lahuta, "FastNS_");
  py::class_<DistanceComputation> DistComp_  (_lahuta, "DistanceComputation_");

  py::class_<TopologyBuildingOptions> Tops_   (_lahuta, "TopologyBuildingOptions");
  py::enum_<ContactComputerType>      CcompT_ (_lahuta, "ContactComputerType");

  DistComp_
    .def_static("distance", py::overload_cast<const Vector<float>&,  const Vector<float>&> (&DistanceComputation::distance<float>))
    .def_static("distance", py::overload_cast<const Vector<double>&, const Vector<double>&>(&DistanceComputation::distance<double>))

    .def_static("distance", py::overload_cast<const std::vector<std::vector<float>> &> (&distance<float>))
    .def_static("distance", py::overload_cast<const std::vector<std::vector<double>> &>(&distance<double>))

    .def_static("distance", py::overload_cast<const std::vector<std::vector<float>> &,  const std::vector<std::vector<float>> &> (&distance<float>))
    .def_static("distance", py::overload_cast<const std::vector<std::vector<double>> &, const std::vector<std::vector<double>> &>(&distance<double>))

    // only double support
    .def_static("search", py::overload_cast<const Matrix<double>&,                        double>( &DistanceComputation::search))
    .def_static("search", py::overload_cast<const Matrix<double>&, const Matrix<double>&, double>( &DistanceComputation::search));


  FastNS_
    .def(py::init<>())
    .def(py::init<const std::vector<RDGeom::Point3D>&,     float>(), py::arg("coords"), py::arg("scale_factor") = 1.1f)
    .def(py::init<const std::vector<std::vector<double>>&, float>(), py::arg("coords"), py::arg("scale_factor") = 1.1f)
    .def(py::init([](py::array_t<double, py::array::c_style | py::array::forcecast> coords_np, float scale_factor) {
      py::buffer_info buf = coords_np.request();

      if (buf.ndim != 2 || buf.shape[1] != 3) throw std::invalid_argument("Input numpy array must have shape (n, 3)");

      size_t n = buf.shape[0];
      double* ptr = static_cast<double*>(buf.ptr);
      std::vector<RDGeom::Point3D> points;
      points.reserve(n);
      for (size_t i = 0; i < n; i++) {
        points.emplace_back(ptr[i * 3], ptr[i * 3 + 1], ptr[i * 3 + 2]);
      }
      return new FastNS(points, scale_factor);

    }))

    .def("build",          &FastNS::build)
    .def("update",         &FastNS::update)
    .def("adaptive_build", &FastNS::adaptive_build)
    .def("self_search",    &FastNS::self_search)
    .def("search",         py::overload_cast<const RDGeom::POINT3D_VECT &>(&FastNS::search, py::const_))
    .def("search",         py::overload_cast<const std::vector<std::vector<double>> &>(&FastNS::search, py::const_))
    .def("get_cutoff",     &FastNS::get_cutoff)

    .def_static("dist_sq", [](const float *a, const float *b) { return FastNS::dist_sq(a, b); })
    .def_static("dist",    [](const float *a, const float *b) { return sqrt(FastNS::dist_sq(a, b)); });

  NSResults_
    .def(py::init<>())
    .def(py::init<const NSResults &>())
    .def(py::init([](py::array_t<int, py::array::c_style | py::array::forcecast> py_pairs, py::array_t<float, py::array::c_style | py::array::forcecast> py_dists) {
      py::buffer_info buf_pairs = py_pairs.request();
      py::buffer_info buf_dists = py_dists.request();

      if (buf_dists.ndim != 1)                            throw std::invalid_argument("distances must be a 1D array");
      if (buf_pairs.ndim != 2 || buf_pairs.shape[1] != 2) throw std::invalid_argument("pairs must be a 2D array with shape (n, 2)");
      if (buf_pairs.shape[0] != buf_dists.shape[0])       throw std::invalid_argument("Number of pairs must match number of distances");

      size_t n = buf_pairs.shape[0];
      Pairs pairs_vec;
      pairs_vec.reserve(n);
      int* pairs_ptr = static_cast<int*>(buf_pairs.ptr);
      for (size_t i = 0; i < n; ++i) {
        pairs_vec.push_back({ pairs_ptr[i * 2], pairs_ptr[i * 2 + 1] });
      }
      float* dists_ptr = static_cast<float*>(buf_dists.ptr);
      std::vector<float> dists_vec(dists_ptr, dists_ptr + n);

      return new NSResults(std::move(pairs_vec), std::move(dists_vec));
    }))
    .def("add",  &NSResults::add)
    .def("size", &NSResults::size)
    .def("add_neighbors", &NSResults::add_neighbors)
    .def("reserve_space", &NSResults::reserve_space)
    .def("filter", py::overload_cast<double>                       (&NSResults::filter, py::const_))
    .def("filter", py::overload_cast<const std::vector<int> &>     (&NSResults::filter, py::const_))
    .def("filter", py::overload_cast<const std::vector<int> &, int>(&NSResults::filter, py::const_))


    .def("clear", &NSResults::clear)
    /*.def("get_pairs",     &NSResults::get_pairs,     py::return_value_policy::reference_internal) */
    /*.def("get_distances", &NSResults::get_distances, py::return_value_policy::reference_internal) */
    .def("get_pairs", [](NSResults &nsresults) {
      auto pairs = nsresults.get_pairs();
      ssize_t n = pairs.size();
      ssize_t dim = 2;
      auto result = py::array_t<int>({n, dim});
      auto buf = result.request();
      int *ptr = static_cast<int *>(buf.ptr);

      for (size_t i = 0; i < n; ++i) {
        ptr[i * 2] = pairs[i].first;
        ptr[i * 2 + 1] = pairs[i].second;
      }

      return result;
    })
    .def("get_distances", [](NSResults &nsresults) {
      auto distances = nsresults.get_distances();
      ssize_t n = distances.size();
      auto result = py::array_t<float>(n);
      auto buf = result.request();
      float *ptr = static_cast<float *>(buf.ptr);

      for (size_t i = 0; i < n; ++i) {
        ptr[i] = distances[i];
      }

      return result;
    })
    .def("__len__", &NSResults::size)
    .def("__iter__", [](const NSResults &nsresults) { return py::make_iterator(nsresults.begin(), nsresults.end()); }, py::keep_alive<0, 1>());


  IntRepr_
      .def(py::init<std::vector<int>, std::vector<int>, std::vector<std::string>,
                    std::vector<int>, std::vector<std::string>,
                    std::vector<std::string>, std::vector<std::vector<float>>>())
      .def_readwrite("atom_indices", &IR::atom_indices)
      .def_readwrite("atom_nums",    &IR::atomic_numbers)
      .def_readwrite("atom_names",   &IR::atom_names)
      .def_readwrite("resids",       &IR::resids)
      .def_readwrite("resnames",     &IR::resnames)
      .def_readwrite("chainlabels",  &IR::chainlabels)
      .def_readwrite("positions",    &IR::positions);


  Residue_
      .def(py::init<>())
      .def(py::init<const std::string &, int, const std::string &, const std::string &>())
      /*.def(py::init<const std::string &, int, const std::string &, bool>())*/
      .def_readwrite("chain_id", &Residue::chain_id)
      .def_readwrite("number",   &Residue::number)
      .def_readwrite("name",     &Residue::name)
      .def_readwrite("alt_loc",  &Residue::alt_loc)
      .def_readwrite("atoms",    &Residue::atoms);


  Residues_
      .def(py::init<const RDKit::RWMol &>(), py::arg("mol"))
      .def_property_readonly("residues", &Residues::get_residues)
      .def("filter", [](const Residues &self, py::function func) {
          // wrap the python callable as a C++ lambda.
          auto pred = [func](const Residue &r) {return func(r).cast<bool>();};
          return self.filter(pred);
      }, py::arg("func"))

      .def("map", [](const Residues &self, py::function func) {
          std::vector<py::object> results;
          for (const auto &r : self.get_residues()) {
              results.push_back(func(r));
          }
          return results;
      }, py::arg("func"))

      .def("__getitem__", [](const Residues &self, size_t i) {if (i >= self.get_residues().size()) throw py::index_error(); return self.get_residues()[i];})
      .def("__iter__",    [](const Residues &self) {return py::make_iterator(self.begin(), self.end());}, py::keep_alive<0, 1>());


  Topology_
      .def_property_readonly("atom_types", &Topology::get_atom_types)
      .def_property_readonly("residues",   [](Topology &top) {return top.get_residues();}, py::return_value_policy::reference)
      .def_property_readonly("rings",      &Topology::get_rings);


  CcompT_
      .value("None",     ContactComputerType::None)
      .value("Arpeggio", ContactComputerType::Arpeggio)
      .value("Molstar",  ContactComputerType::Molstar);


  Tops_
      .def(py::init<>())
      .def_readwrite("atom_typing_method",   &TopologyBuildingOptions::atom_typing_method)
      .def_readwrite("bonded_search_cutoff", &TopologyBuildingOptions::cutoff);


  Luni
    /*.def(py::init<std::string>())*/
      /*.def(py::init<const IR &>())*/ // FIX: this uses the `create` factory function now
      .def(py::init<std::string>(), py::arg("file_name"))
      .def("build_topology", &Luni::build_topology, py::arg("t_opts") = std::nullopt)
      .def("find_neighbors", [](class Luni &luni, double cutoff, int res_dif){
          auto grid = FastNS(luni.get_conformer().getPositions());
          if (!grid.build(cutoff)) {
            spdlog::error("Failed to update the grid with the given cutoff");
            return NSResults();
          }
          auto ns = grid.self_search();
          if (res_dif > 0) {
            ns = luni.remove_adjascent_residueid_pairs(ns, res_dif);
          }
          return ns;
      })
      .def("filter",           &Luni::filter)
      .def("parse_expression", &Luni::parse_expression)
      .def("get_atom_types",   &Luni::get_atom_types)
      .def("get_rings",        &Luni::get_rings)
      .def("filter_luni",      &Luni::filter_luni)
      /*.def("match_smarts_string", &Luni::match_smarts_string)*/

      .def("get_positions", [](class Luni &luni, int id = -1)    {auto values = luni.get_conformer(id).getPositions(); return coordinates(values);})
      .def_property_readonly("positions",   [](class Luni &luni) {auto values = luni.get_conformer()  .getPositions(); return coordinates(values);})

      .def_property_readonly("indices",     [](class Luni &luni) {auto values = luni.indices();        return int_array(values);})
      .def_property_readonly("atom_nums",   [](class Luni &luni) {auto values = luni.atomic_numbers(); return int_array(values);})
      .def_property_readonly("resids",      [](class Luni &luni) {auto values = luni.resids();         return int_array(values);})
      .def_property_readonly("resindices",  [](class Luni &luni) {auto values = luni.resindices();     return int_array(values);})
      .def_property_readonly("names",       [](class Luni &luni) {auto values = luni.names();          return string_array(values);})
      .def_property_readonly("symbols",     [](class Luni &luni) {auto values = luni.symbols();        return string_array(values);})
      .def_property_readonly("elements",    [](class Luni &luni) {auto values = luni.elements();       return string_array(values);})
      .def_property_readonly("resnames",    [](class Luni &luni) {auto values = luni.resnames();       return string_array(values);})
      .def_property_readonly("chainlabels", [](class Luni &luni) {auto values = luni.chainlabels();    return string_array(values);})

      .def_property_readonly("n_atoms",     [](class Luni &luni) {return luni.n_atoms();})
      .def_property_readonly("file_name",   [](class Luni &luni) {return luni.file_name_.c_str();})

       .def("get_topology", [](class Luni &luni) -> const Topology& { return luni.get_topology(); }, py::return_value_policy::reference)
       /*.def("at", [](class Luni &luni) -> auto { auto &v = luni.get_topology(); return v.atom_types; })*/
      .def("get_atom", &Luni::get_atom, py::return_value_policy::reference)

      .def("count_unique", py::overload_cast<const std::vector<int> &>        (&Luni::count_unique))
      .def("count_unique", py::overload_cast<const std::vector<std::string> &>(&Luni::count_unique))
      .def("find_elements", &Luni::find_elements)
      .def("factorize",     &Luni::factorize)

      // FIX: these will be moved to the topology class
      .def("assign_molstar_atom_types",  &Luni::assign_molstar_atom_types)
      .def("assign_arpeggio_atom_types", &Luni::assign_arpeggio_atom_types);

  // FIX: We provide a NumPy-based interface, so we should, perhaps, not expose the RDKit Poin3D class.
  py::class_<RDGeom::Point3D>(_lahuta, "Point3D")
      .def(py::init<>())
      .def(py::init<double, double, double>())
      .def(py::init([](py::array_t<double, py::array::c_style | py::array::forcecast> coords) {
        py::buffer_info buf = coords.request();
        if (buf.ndim != 1 || buf.shape[0] != 3)
          throw std::invalid_argument("Input numpy array must have shape (3)");
        double* ptr = static_cast<double*>(buf.ptr);
        return RDGeom::Point3D(ptr[0], ptr[1], ptr[2]);
      }))
      .def_readwrite("x", &RDGeom::Point3D::x)
      .def_readwrite("y", &RDGeom::Point3D::y)
      .def_readwrite("z", &RDGeom::Point3D::z)
      .def("__str__",  [](const RDGeom::Point3D &p) { return "(" + std::to_string(p.x) + ", " + std::to_string(p.y) + ", " + std::to_string(p.z) + ")"; });


  Atom_
    .def_property_readonly("atomic_number",    &RDKit::Atom::getAtomicNum)
    .def_property_readonly("symbol",           &RDKit::Atom::getSymbol)
    .def_property_readonly("idx",              &RDKit::Atom::getIdx)
    .def_property_readonly("degree",           &RDKit::Atom::getDegree)
    .def_property_readonly("formal_charge",    &RDKit::Atom::getFormalCharge)
    .def_property_readonly("hybridization",    &RDKit::Atom::getHybridization)
    .def_property_readonly("is_aromatic",      &RDKit::Atom::getIsAromatic)
    .def_property_readonly("mass",             &RDKit::Atom::getMass)
    .def_property_readonly("num_explicit_Hs",  &RDKit::Atom::getNumExplicitHs)
    .def_property_readonly("num_implicit_Hs",  &RDKit::Atom::getNumImplicitHs)
    .def_property_readonly("total_num_Hs",     &RDKit::Atom::getTotalNumHs)
    .def_property_readonly("total_valence",    &RDKit::Atom::getTotalValence)
    .def_property_readonly("explicit_valence", &RDKit::Atom::getExplicitValence)
    .def_property_readonly("implicit_valence", &RDKit::Atom::getImplicitValence);

}

PYBIND11_MODULE(_lahuta, m) {
  m.doc() = "lahuta: A Python binding for the Lahuta library";

  py::class_<RDKit::RWMol> lahutaRWMol(m, "RWMol");
  py::class_<RDKit::Conformer> lahutaConformer(m, "Conformer");

  bind_atom_types(m);
  bind_contacts(m);
  bind(m);
  bind_common(m);
  bind_entities(m);
}
