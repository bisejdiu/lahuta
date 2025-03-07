#include <pybind11/complex.h> // FIX: is this needed?
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/RDKitBase.h>

#include "at.hpp"
#include "atom_types.hpp"
#include "common.hpp"
#include "contacts.hpp"
#include "lahuta.hpp"
#include "entities.hpp"
#include "neighbors.hpp"
#include "nsgrid.hpp"
#include "residues.hpp"
#include "topology.hpp"
/*#include "nn.hpp"*/
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
  py::class_<Luni>                Luni       (_lahuta, "LahutaCPP");
  py::class_<Topology>            Topology_  (_lahuta, "Topology_");
  py::class_<Residue>             Residue_   (_lahuta, "Residue");
  py::class_<Residues>            Residues_  (_lahuta, "Residues_");
  py::class_<IR>                  IntRepr_   (_lahuta, "IR");
  py::class_<NSResults>           NSResults_ (_lahuta, "NSResults_");
  py::class_<FastNS>              FastNS_    (_lahuta, "FastNS_");
  py::class_<DistanceComputation> DistComp_  (_lahuta, "DistanceComputation_");


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
    .def(py::init<const std::vector<RDGeom::Point3D>&, float>(), py::arg("coords"), py::arg("scale_factor") = 1.1f)
    .def(py::init<const std::vector<std::vector<double>>&, float>(), py::arg("coords"), py::arg("scale_factor") = 1.1f)
    .def(py::init([](py::array_t<double, py::array::c_style | py::array::forcecast> coords_np, float scale_factor) {
      py::buffer_info buf = coords_np.request();
      if (buf.ndim != 2 || buf.shape[1] != 3)
        throw std::invalid_argument("Input numpy array must have shape (n, 3)");
      size_t n = buf.shape[0];
      double* ptr = static_cast<double*>(buf.ptr);
      std::vector<RDGeom::Point3D> points;
      points.reserve(n);
      for (size_t i = 0; i < n; i++) {
        points.emplace_back(ptr[i * 3], ptr[i * 3 + 1], ptr[i * 3 + 2]);
      }
      return new FastNS(points, scale_factor);
    }), py::arg("coords"), py::arg("scale_factor") = 1.1f)

    .def("build",          &FastNS::build)
    .def("update",         &FastNS::update)
    .def("adaptive_build", &FastNS::adaptive_build)
    .def("self_search",    &FastNS::self_search)
    /*.def("search",         &FastNS::search)*/
    .def("search",       py::overload_cast<const RDGeom::POINT3D_VECT &>(&FastNS::search, py::const_))
    .def("search",       py::overload_cast<const std::vector<std::vector<double>> &>(&FastNS::search, py::const_))
    .def("get_cutoff",     &FastNS::get_cutoff)

    // template <typename T>
    // static inline T dist_sq(const T* __restrict a, const T* __restrict b) {
    //     T dx = a[0] - b[0];
    //     T dy = a[1] - b[1];
    //     T dz = a[2] - b[2];
    //     return dx * dx + dy * dy + dz * dz;
    // }

    /*.def_static("dist_sq", &FastNS::dist_sq)*/
    .def_static("dist_sq", [](const float *a, const float *b) { return FastNS::dist_sq(a, b); })
    .def_static("dist", [](const float *a, const float *b) { return sqrt(FastNS::dist_sq(a, b)); });

  NSResults_
    // Default constructor
    .def(py::init<>())
    // Expose copy constructor
    .def(py::init<const NSResults &>())
    // Custom constructor: construct NSResults from two Python lists.
    .def(py::init([](py::array_t<int, py::array::c_style | py::array::forcecast> py_pairs, py::array_t<float, py::array::c_style | py::array::forcecast> py_dists) {
      // Request buffer info
      py::buffer_info buf_pairs = py_pairs.request();
      py::buffer_info buf_dists = py_dists.request();

      // Check that pairs is 2D with second dimension 2
      if (buf_pairs.ndim != 2 || buf_pairs.shape[1] != 2)
        throw std::invalid_argument("pairs must be a 2D array with shape (n, 2)");
      // Check that distances is 1D and length matches the number of pairs
      if (buf_dists.ndim != 1)
        throw std::invalid_argument("distances must be a 1D array");
      if (buf_pairs.shape[0] != buf_dists.shape[0])
        throw std::invalid_argument("Number of pairs must match number of distances");

      size_t n = buf_pairs.shape[0];
      Pairs pairs_vec;
      pairs_vec.reserve(n);
      int* pairs_ptr = static_cast<int*>(buf_pairs.ptr);
      for (size_t i = 0; i < n; ++i) {
        // Each row in pairs is two consecutive ints
        pairs_vec.push_back({ pairs_ptr[i * 2], pairs_ptr[i * 2 + 1] });
      }
      float* dists_ptr = static_cast<float*>(buf_dists.ptr);
      std::vector<float> dists_vec(dists_ptr, dists_ptr + n);

      return new NSResults(std::move(pairs_vec), std::move(dists_vec));
    }))
    .def("add", &NSResults::add)
    .def("add_neighbors", &NSResults::add_neighbors)
    .def("reserve_space", &NSResults::reserve_space)
    .def("size", &NSResults::size)
    .def("__len__", &NSResults::size)  // Allow use of len() in Python.
    // Bind the overloads for filter.
    /*.def("filter", (NSResults (NSResults::*)(double) const) &NSResults::filter)*/
    /*.def("filter", (NSResults (NSResults::*)(const std::vector<int> &) const) &NSResults::filter)*/

    .def("filter", py::overload_cast<double>                  (&NSResults::filter, py::const_))
    .def("filter", py::overload_cast<const std::vector<int> &>(&NSResults::filter, py::const_))

    .def("clear", &NSResults::clear)
    /*.def("get_pairs", &NSResults::get_pairs, py::return_value_policy::reference_internal) // this returns list of pairs, but we should return a numpy array*/
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
    /*.def("get_distances", &NSResults::get_distances, py::return_value_policy::reference_internal) // this returns list of distances, but we should return a numpy array*/
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
    // Enable iteration support.
    .def("__iter__", [](const NSResults &nsresults) {
      return py::make_iterator(nsresults.begin(), nsresults.end());
    }, py::keep_alive<0, 1>());

  ////!!!! Recs_
  ////!!!!     .def(py::init<>())
  ////!!!!     /*.def_readwrite("rings", &RingDataVec::rings)*/

  ////!!!!     .def_property_readonly("rings", [](RingEntityCollection &rdv) { return rdv.get_data(); })
  ////!!!!     .def("compute_angles", &RingEntityCollection::compute_angles)
  ////!!!!     /*.def("root_atom_ids", &RingDataVec::root_atom_ids)*/
  ////!!!!     /*.def("compute_angles",*/
  ////!!!!     /*     [](RingDataVec &rdv, const std::vector<std::vector<double>>
  ////!!!!        &points) {*/
  ////!!!!     /*       return rdv.compute_angles(points);*/
  ////!!!!     /*     })*/
  ////!!!!     // .def_property_readonly("centers", &RingDataVec::centers)
  ////!!!!     // .def_property_readonly("centers", [](RingDataVec &rdv) {
  ////!!!!     //   auto centers = rdv.centers();
  ////!!!!     //   ssize_t n = centers.size();
  ////!!!!     //   ssize_t dim = 3;
  ////!!!!     //   auto result = py::array_t<double>({n, dim});
  ////!!!!     //   auto buf = result.request();
  ////!!!!     //   double *ptr = static_cast<double *>(buf.ptr);
  ////!!!!     //
  ////!!!!     //   for (size_t i = 0; i < n; ++i) {
  ////!!!!     //     for (size_t j = 0; j < dim; ++j) {
  ////!!!!     //       ptr[i * dim + j] = centers[i][j];
  ////!!!!     //     }
  ////!!!!     //   }
  ////!!!!     //
  ////!!!!     //   return result;
  ////!!!!     // })
  ////!!!!     .def_property_readonly("centers",
  ////!!!!                            [](RingEntityCollection &rdv) {
  ////!!!!                              ssize_t n = rdv.data.size();
  ////!!!!                              ssize_t dim = 3;
  ////!!!!                              auto result = py::array_t<double>({n, dim});
  ////!!!!                              auto buf = result.request();
  ////!!!!                              double *ptr = static_cast<double *>(buf.ptr);

  ////!!!!                              for (size_t i = 0; i < n; ++i) {
  ////!!!!                                ptr[i * dim] = rdv.data[i].center.x;
  ////!!!!                                ptr[i * dim + 1] = rdv.data[i].center.y;
  ////!!!!                                ptr[i * dim + 2] = rdv.data[i].center.z;
  ////!!!!                              }

  ////!!!!                              return result;
  ////!!!!                            })
  ////!!!!     .def_property_readonly("norm",
  ////!!!!                            [](RingEntityCollection &rdv) {
  ////!!!!                              ssize_t n = rdv.data.size();
  ////!!!!                              ssize_t dim = 3;
  ////!!!!                              auto result = py::array_t<double>({n, dim});
  ////!!!!                              auto buf = result.request();
  ////!!!!                              double *ptr = static_cast<double *>(buf.ptr);

  ////!!!!                              for (size_t i = 0; i < n; ++i) {
  ////!!!!                                ptr[i * dim] = rdv.data[i].norm.x;
  ////!!!!                                ptr[i * dim + 1] = rdv.data[i].norm.y;
  ////!!!!                                ptr[i * dim + 2] = rdv.data[i].norm.z;
  ////!!!!                              }

  ////!!!!                              return result;
  ////!!!!                            });

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
      .def(py::init([](class Luni &luni, const Pairs &&pairs, const Distances &&dists, bool is_sorted) {
        return Neighbors<AtomAtomPair>(luni, std::move(pairs), std::move(dists), is_sorted);
      }))
      .def_property_readonly("data", [](Neighbors<AtomAtomPair> &self) { return self.get_data(); })
      /*.def("data", &Neighbors<AtomAtomPair>::data) */
      .def("size", &Neighbors<AtomAtomPair>::size)
      // .def("type_filter", &Neighbors<AtomAtomPair>::type_filter)
      .def("type_filter", [](Neighbors<AtomAtomPair> &self, AtomType type, int partner) {return self.type_filter(type, partner);})
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
      // FIX: we can expose the two types of declarations but it's better to handle the calls on the Python side
      .def("intersection_x", [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) {return Neighbors<AtomAtomPair>::intersection(self.get_data(), other.get_data());})
      .def("intersection",         [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) {return self.intersection(other);})
      .def("difference",           [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) { return self.difference(other);})
      .def("union",                [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) { return self.union_(other);})
      .def("symmetric_difference", [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) { return self.symmetric_difference(other); })

      .def("__add__", [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) { return self.union_(other); })
      .def("__sub__", [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) { return self.difference(other); })
      .def("__xor__", [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) { return self.symmetric_difference(other);})
      .def("__and__", [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) { return self.intersection(other); })
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
        return Neighbors<AtomRingPair>(luni, std::move(pairs), std::move(dists), is_sorted);
      }))
      // .def(py::init([](class Luni &luni, Pairs pairs, Distances dists,
      //                  ContactType c,
      //                  bool is_sorted = false) {
      //   return Neighbors<AtomAtomPair>(luni, pairs, dists, c, is_sorted);
      // }))

      .def_property_readonly("data", [](Neighbors<AtomRingPair> &self) { return self.get_data(); })
      .def("size", &Neighbors<AtomRingPair>::size)
      .def("type_filter", [](Neighbors<AtomRingPair> &self, AtomType type, int partner) { return self.type_filter(type, partner);})
      // .def("type_filter", &Neighbors<AtomRingPair>::type_filter)

      .def("get_pairs",
           [](Neighbors<AtomRingPair> &nsr) {
             const auto &pairs = nsr.get_pairs();
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

      // FIX: here we do not return numpy arrays
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

      .def("__add__", [](Neighbors<AtomRingPair> &self, Neighbors<AtomRingPair> &other) {return self.union_(other);})
      .def("__sub__", [](Neighbors<AtomRingPair> &self, Neighbors<AtomRingPair> &other) {return self.difference(other);})
      .def("__and__", [](Neighbors<AtomRingPair> &self, Neighbors<AtomRingPair> &other) {return self.intersection(other);})
      .def("__xor__", [](Neighbors<AtomRingPair> &self, Neighbors<AtomRingPair> &other) {return self.symmetric_difference(other);})
      .def("get_luni",  &Neighbors<AtomRingPair>::get_luni);

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
      .def_property_readonly("residues",   [](Topology &top) {return top.residues.get();}, py::return_value_policy::reference)
      .def_property_readonly("rings",      &Topology::get_rings);


  Luni
    /*.def(py::init<std::string>())*/
      .def(py::init<const IR &>())
      .def(py::init<std::string, std::optional<ContactComputerType>>(), py::arg("file_name"), py::arg("c_type") = std::nullopt)
      .def(py::init([](std::string file_name, int c) { return new class Luni(file_name, static_cast<ContactComputerType>(c)); }), py::arg("file_name"), py::arg("contact_type"))
      .def("find_neighbors", &Luni::find_neighbors)
      /*.def("find_neighbors_aa", &Luni::find_neighbors<AtomAtomPair>)*/
      /*.def("find_neighbors_ar", &Luni::find_neighbors<AtomRingPair>)*/
      /*.def("fn_aa", &Luni::find_neighbors<AtomAtomPair>)*/
      /*.def("fn_ar", &Luni::find_neighbors<AtomRingPair>)*/
      /*.def("find_ring_neighbors", &Luni::find_ring_neighbors)*/
      /*.def("find_neighbors",*/
      /*     [](class Luni &luni, float cutoff) {*/
      /*       return luni.find_neighbors_opt(cutoff);*/
      /*     })*/

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
      .def_property_readonly("file_name",   [](class Luni &luni) {return luni.file_name.c_str();})

       .def("get_topology", [](class Luni &luni) -> const Topology& { return luni.get_topology(); }, py::return_value_policy::reference)
       /*.def("at", [](class Luni &luni) -> auto { auto &v = luni.get_topology(); return v.atom_types; })*/

      // FIX: either provide the property method or the function, but not both. I think the former are better.
      .def("get_indices",        &Luni::indices)
      .def("get_atomic_numbers", &Luni::atomic_numbers)
      .def("get_names",          &Luni::names)
      .def("get_symbols",        &Luni::symbols)
      .def("get_elements",       &Luni::elements)
      .def("get_resnames",       &Luni::resnames)
      .def("get_resids",         &Luni::resids)
      .def("get_resindices",     &Luni::resindices)
      .def("get_chainlabels",    &Luni::chainlabels)

      .def("get_atom", &Luni::get_atom, py::return_value_policy::reference)

      .def("count_unique", py::overload_cast<const std::vector<int> &>        (&Luni::count_unique))
      .def("count_unique", py::overload_cast<const std::vector<std::string> &>(&Luni::count_unique))
      .def("find_elements", &Luni::find_elements)
      .def("factorize",     &Luni::factorize)

       // FIX: selection should be done using an enum
      .def("assign_molstar_atom_types",  &Luni::assign_molstar_atom_types)
      .def("assign_arpeggio_atom_types", &Luni::assign_arpeggio_atom_types)

      .def("get_cutoff", &Luni::get_cutoff);

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

  py::class_<RDKit::Atom> Atom_(_lahuta, "Atom");

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
