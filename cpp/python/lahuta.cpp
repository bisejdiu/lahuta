#include <pybind11/complex.h> // FIX: is this needed?
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/RDKitBase.h>

#include "common.hpp"
#include "contacts.hpp"
#include "at.hpp"
#include "lahuta.hpp"
#include "neighbors.hpp"
#include "nsgrid.hpp"
/*#include "nn.hpp"*/

namespace py = pybind11;
using namespace lahuta;

void bind(py::module &_lahuta) {
  py::class_<Luni> Luni(_lahuta, "LahutaCPP");

  py::class_<RingData>(_lahuta, "RingData")
      .def(py::init<>())
      /*.def_readwrite("atom_ids", &RingData::atom_ids)*/
      .def("compute_angle", &RingData::compute_angle)
      /*.def("compute_angle2", &RingData::compute_angle2)*/
      .def_property_readonly(
          "center", [](RingData &rd) { return point3d_to_pyarray(rd.center); })
      .def_property_readonly(
          "norm1", [](RingData &rd) { return point3d_to_pyarray(rd.norm); });

  py::class_<RingDataVec>(_lahuta, "RingDataVec")
      .def(py::init<>())
      .def_readwrite("rings", &RingDataVec::rings)
      .def("compute_angles", &RingDataVec::compute_angles)
      .def("root_atom_ids", &RingDataVec::root_atom_ids)
      /*.def("compute_angles",*/
      /*     [](RingDataVec &rdv, const std::vector<std::vector<double>>
         &points) {*/
      /*       return rdv.compute_angles(points);*/
      /*     })*/
      // .def_property_readonly("centers", &RingDataVec::centers)
      // .def_property_readonly("centers", [](RingDataVec &rdv) {
      //   auto centers = rdv.centers();
      //   ssize_t n = centers.size();
      //   ssize_t dim = 3;
      //   auto result = py::array_t<double>({n, dim});
      //   auto buf = result.request();
      //   double *ptr = static_cast<double *>(buf.ptr);
      //
      //   for (size_t i = 0; i < n; ++i) {
      //     for (size_t j = 0; j < dim; ++j) {
      //       ptr[i * dim + j] = centers[i][j];
      //     }
      //   }
      //
      //   return result;
      // })
      .def_property_readonly("centers",
                             [](RingDataVec &rdv) {
                               ssize_t n = rdv.rings.size();
                               ssize_t dim = 3;
                               auto result = py::array_t<double>({n, dim});
                               auto buf = result.request();
                               double *ptr = static_cast<double *>(buf.ptr);

                               for (size_t i = 0; i < n; ++i) {
                                 ptr[i * dim] = rdv.rings[i].center.x;
                                 ptr[i * dim + 1] = rdv.rings[i].center.y;
                                 ptr[i * dim + 2] = rdv.rings[i].center.z;
                               }

                               return result;
                             })
      .def_property_readonly("norm",
                             [](RingDataVec &rdv) {
                               ssize_t n = rdv.rings.size();
                               ssize_t dim = 3;
                               auto result = py::array_t<double>({n, dim});
                               auto buf = result.request();
                               double *ptr = static_cast<double *>(buf.ptr);

                               for (size_t i = 0; i < n; ++i) {
                                 ptr[i * dim] = rdv.rings[i].norm.x;
                                 ptr[i * dim + 1] = rdv.rings[i].norm.y;
                                 ptr[i * dim + 2] = rdv.rings[i].norm.z;
                               }

                               return result;
                             });

  // py::enum_<ContactType>(_lahuta, "ContactType")
  //     .value("AtomAtom", ContactType::AtomAtom)
  //     .value("AtomRing", ContactType::AtomRing);

  py::class_<AtomAtomPair>(_lahuta, "AtomAtomPair")
      .def(py::init<int, int, float>())
      .def_readonly("i", &AtomAtomPair::i)
      .def_readonly("j", &AtomAtomPair::j)
      .def_readonly("d", &AtomAtomPair::d)
      .def("get_pair", &AtomAtomPair::get_pair)
      /*.def("get_i", &AtomAtomPair::get_i)*/
      /*.def("get_j", &AtomAtomPair::get_j)*/
      .def("names", &AtomAtomPair::names);

  py::class_<AtomRingPair>(_lahuta, "AtomRingPair")
      .def(py::init<int, int, float>())
      .def_readonly("i", &AtomRingPair::i)
      .def_readonly("j", &AtomRingPair::j)
      .def_readonly("d", &AtomRingPair::d)
      .def("get_pair", &AtomRingPair::get_pair)
      /*.def("get_i", &AtomRingPair::get_i)*/
      /*.def("get_j", &AtomRingPair::get_j)*/
      .def("names", &AtomRingPair::names);

  py::class_<Neighbors<AtomAtomPair>>(_lahuta, "AtomAtomNeighbors")
      .def(py::init(
          [](class Luni &luni, std::vector<AtomAtomPair> data, bool is_sorted) {
            return Neighbors<AtomAtomPair>(luni, data, is_sorted);
          }))
      .def(py::init([](class Luni &luni, const Pairs &&pairs,
                       const Distances &&dists, bool is_sorted) {
        return Neighbors<AtomAtomPair>(luni, std::move(pairs), std::move(dists),
                                       is_sorted);
      }))
      .def_property_readonly(
          "data", [](Neighbors<AtomAtomPair> &self) { return self.get_data(); })
      /*.def("data", &Neighbors<AtomAtomPair>::data) */
      .def("size", &Neighbors<AtomAtomPair>::size)
      // .def("type_filter", &Neighbors<AtomAtomPair>::type_filter)
      .def("type_filter",
           [](Neighbors<AtomAtomPair> &self, AtomType type, int partner) {
             return self.type_filter(type, partner);
           })
      .def("get_pairs",
           [](Neighbors<AtomAtomPair> &nsr) {
             auto pairs = nsr.get_pairs();
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
      .def("get_distances_sq", &Neighbors<AtomAtomPair>::get_distances)
      .def("get_distances",
           [](Neighbors<AtomAtomPair> &nsr) {
             auto distances = nsr.get_distances();
             ssize_t n = distances.size();
             auto result = py::array_t<float>(n);
             auto buf = result.request();
             float *ptr = static_cast<float *>(buf.ptr);

             for (size_t i = 0; i < n; ++i) {
               // ptr[i] = distances[i];
               ptr[i] = sqrt(distances[i]);
             }

             return result;
           })
      // FIX: we can expose the two types of declarations but it's better to
      // handle the calls on the Python side
      .def("intersection_x",
           [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) {
             return Neighbors<AtomAtomPair>::intersection(self.get_data(),
                                                          other.get_data());
           })
      .def("intersection",
           [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) {
             return self.intersection(other);
           })
      .def("difference",
           [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) {
             return self.difference(other);
           })
      .def("union",
           [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) {
             return self.union_(other);
           })
      .def("symmetric_difference",
           [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) {
             return self.symmetric_difference(other);
           })
      .def("__add__",
           [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) {
             return self.union_(other);
           })
      .def("__sub__",
           [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) {
             return self.difference(other);
           })
      .def("__xor__",
           [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) {
             return self.symmetric_difference(other);
           })
      .def("__and__",
           [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) {
             return self.intersection(other);
           })
      .def("get_luni", &Neighbors<AtomAtomPair>::get_luni);

  py::class_<Neighbors<AtomRingPair>>(_lahuta, "AtomRingNeighbors")
      .def(py::init([](class Luni &luni, std::vector<AtomRingPair> data,
                       /*const AtomRingPair::RefType &ctx,*/
                       bool is_sorted = false) {
        return Neighbors<AtomRingPair>(luni, data, is_sorted);
      }))
      .def(py::init([](class Luni &luni, Pairs pairs, Distances dists,
                       /*const AtomRingPair::RefType &ctx,*/
                       bool is_sorted = false) {
        /*return Neighbors<AtomRingPair>(luni, pairs, dists, is_sorted);*/
        return Neighbors<AtomRingPair>(luni, std::move(pairs), std::move(dists),
                                       is_sorted);
      }))
      // .def(py::init([](class Luni &luni, Pairs pairs, Distances dists,
      //                  ContactType c,
      //                  bool is_sorted = false) {
      //   return Neighbors<AtomAtomPair>(luni, pairs, dists, c, is_sorted);
      // }))
      .def_property_readonly(
          "data", [](Neighbors<AtomRingPair> &self) { return self.get_data(); })
      .def("size", &Neighbors<AtomRingPair>::size)
      // .def("type_filter", &Neighbors<AtomRingPair>::type_filter)
      .def("type_filter",
           [](Neighbors<AtomRingPair> &self, AtomType type, int partner) {
             return self.type_filter(type, partner);
           })

      .def("get_pairs",
           [](Neighbors<AtomRingPair> &nsr) {
             auto pairs = nsr.get_pairs();
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
      .def("get_distances_sq", &Neighbors<AtomRingPair>::get_distances)
      .def("get_distances",
           [](Neighbors<AtomRingPair> &nsr) {
             auto distances = nsr.get_distances();
             ssize_t n = distances.size();
             auto result = py::array_t<float>(n);
             auto buf = result.request();
             float *ptr = static_cast<float *>(buf.ptr);

             for (size_t i = 0; i < n; ++i) {
               // ptr[i] = distances[i];
               ptr[i] = sqrt(distances[i]);
             }

             return result;
           })
      .def("__add__",
           [](Neighbors<AtomRingPair> &self, Neighbors<AtomRingPair> &other) {
             return self.union_(other);
           })
      .def("__sub__",
           [](Neighbors<AtomRingPair> &self, Neighbors<AtomRingPair> &other) {
             return self.difference(other);
           })
      .def("__xor__",
           [](Neighbors<AtomRingPair> &self, Neighbors<AtomRingPair> &other) {
             return self.symmetric_difference(other);
           })
      .def("__and__",
           [](Neighbors<AtomRingPair> &self, Neighbors<AtomRingPair> &other) {
             return self.intersection(other);
           })
      .def("get_luni", &Neighbors<AtomRingPair>::get_luni);

  py::class_<IR>(_lahuta, "IR")
      .def(
          py::init<std::vector<int>, std::vector<int>, std::vector<std::string>,
                   std::vector<int>, std::vector<std::string>,
                   std::vector<std::string>, std::vector<std::vector<float>>>())
      .def_readwrite("atom_indices", &IR::atom_indices)
      .def_readwrite("atomic_numbers", &IR::atomic_numbers)
      .def_readwrite("atom_names", &IR::atom_names)
      .def_readwrite("resids", &IR::resids)
      .def_readwrite("resnames", &IR::resnames)
      .def_readwrite("chainlabels", &IR::chainlabels)
      .def_readwrite("positions", &IR::positions);

  Luni.def(py::init<std::string>())
      .def(py::init<const IR &>())
      .def_property_readonly(
          "file_name", [](class Luni &luni) { return luni.file_name.c_str(); })
      .def("find_neighbors", &Luni::find_neighbors)
      /*.def("find_neighbors_aa", &Luni::find_neighbors<AtomAtomPair>)*/
      /*.def("find_neighbors_ar", &Luni::find_neighbors<AtomRingPair>)*/
      /*.def("fn_aa", &Luni::find_neighbors<AtomAtomPair>)*/
      /*.def("fn_ar", &Luni::find_neighbors<AtomRingPair>)*/
      .def("find_ring_neighbors", &Luni::find_ring_neighbors)
      /*.def("find_neighbors",*/
      /*     [](class Luni &luni, float cutoff) {*/
      /*       return luni.find_neighbors_opt(cutoff);*/
      /*     })*/
      .def("filter", &Luni::filter)
      .def("parse_expression", &Luni::parse_expression)
      .def("get_atom_types", &Luni::get_atom_types)
      .def("get_rings", &Luni::get_rings)
      /*.def("match_smarts_string", &Luni::match_smarts_string)*/
      .def("filter_luni", &Luni::filter_luni)
      /*.def("get_n_atoms", &Luni::n_atoms)*/
      .def_property_readonly("n_atoms",
                             [](class Luni &luni) { return luni.n_atoms(); })

      .def("get_positions", // need to return list[list[float]]
           [](class Luni &luni) {
             std::vector<std::vector<double>> positions;
             auto coords = luni.get_molecule().getConformer().getPositions();
             for (const auto &coord : coords) {
               positions.push_back({coord.x, coord.y, coord.z});
             }
             return positions;
           })
      .def_property_readonly(
          "positions",
          [](class Luni &luni) {
            return coordinates(
                luni.get_molecule().getConformer().getPositions());
          })

      .def("get_indices", &Luni::indices)
      .def_property_readonly("indices",
                             [](class Luni &luni) {
                               auto indices = luni.indices();
                               return int_array(indices);
                             })
      .def("get_atomic_numbers", &Luni::atomic_numbers)
      .def_property_readonly("atomic_numbers",
                             [](class Luni &luni) {
                               auto atomic_numbers = luni.atomic_numbers();
                               return int_array(atomic_numbers);
                             })
      .def("get_names", &Luni::names)
      .def_property_readonly("names",
                             [](class Luni &luni) {
                               auto names = luni.names();
                               return string_array(names);
                             })
      .def("get_symbols", &Luni::symbols)
      .def_property_readonly("symbols",
                             [](class Luni &luni) {
                               auto symbols = luni.symbols();
                               return string_array(symbols);
                             })
      .def("get_elements", &Luni::elements)
      .def_property_readonly("elements",
                             [](class Luni &luni) {
                               auto elements = luni.elements();
                               return string_array(elements);
                             })
      .def("get_resnames", &Luni::resnames)
      .def_property_readonly("resnames",
                             [](class Luni &luni) {
                               auto resnames = luni.resnames();
                               return string_array(resnames);
                             })
      .def("get_resids", &Luni::resids)
      .def_property_readonly("resids",
                             [](class Luni &luni) {
                               auto resids = luni.resids();
                               return int_array(resids);
                             })
      .def("get_resindices", &Luni::resindices)
      .def_property_readonly("resindices",
                             [](class Luni &luni) {
                               auto resindices = luni.resindices();
                               return int_array(resindices);
                             })
      .def("get_chainlabels", &Luni::chainlabels)
      .def_property_readonly("chainlabels",
                             [](class Luni &luni) {
                               auto chainlabels = luni.chainlabels();
                               return string_array(chainlabels);
                             })

      .def("count_unique",
           py::overload_cast<const std::vector<int> &>(&Luni::count_unique))
      .def("count_unique", py::overload_cast<const std::vector<std::string> &>(
                               &Luni::count_unique))
      .def("find_elements", &Luni::find_elements)
      .def("factorize", &Luni::factorize)
      
      .def("test_find_neighbors", &Luni::test_find_neighbors)
      .def("assign_molstar_atom_types", &Luni::assign_molstar_atom_types)
      .def("assign_arpeggio_atom_types", &Luni::assign_arpeggio_atom_types)

      .def("get_cutoff", &Luni::get_cutoff);

  py::class_<RDGeom::Point3D>(_lahuta, "Point3D")
      .def(py::init<>())
      .def_readwrite("x", &RDGeom::Point3D::x)
      .def_readwrite("y", &RDGeom::Point3D::y)
      .def_readwrite("z", &RDGeom::Point3D::z);
}

PYBIND11_MODULE(_lahuta, m) {
  m.doc() = "lahuta: A Python binding for the Lahuta library";

  py::class_<RDKit::RWMol> lahutaRWMol(m, "RWMol");
  py::class_<RDKit::Conformer> lahutaConformer(m, "Conformer");

  bind_atom_types(m);
  bind_contacts(m);
  bind(m);
  bind_common(m);
}
