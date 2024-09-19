#include "lahuta.hpp"
#include "neighbors.hpp"
#include "nsgrid.hpp"
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/RDKitBase.h>

#include "at.hpp"

namespace py = pybind11;
using namespace lahuta;

py::array_t<float> coordinates(const RDGeom::POINT3D_VECT &coords) {
  if (coords.empty() || coords[0].dimension() != 3) {
    throw std::runtime_error(
        "Invalid input: expected non-empty vector of 3D coordinates");
  }

  size_t n_atoms = coords.size();
  size_t n_dims = coords[0].dimension();

  auto result = py::array_t<float>({n_atoms, n_dims});
  auto buf = result.request();
  float *ptr = static_cast<float *>(buf.ptr);

  for (size_t i = 0; i < n_atoms; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      ptr[i * 3 + j] = coords[i][j];
    }
  }

  return result;
}

auto point3d_to_pyarray = [](const RDGeom::Point3D &vec) {
  auto result = py::array_t<float>(3);
  auto buf = result.request();
  float *ptr = static_cast<float *>(buf.ptr);
  ptr[0] = vec.x;
  ptr[1] = vec.y;
  ptr[2] = vec.z;
  return result;
};

void bind(py::module &_lahuta) {
  py::class_<Luni> Luni(_lahuta, "Luni");

  // .def("get_molecule",
  //      (RDKit::RWMol & (GemmiSource::*)()) & GemmiSource::get_molecule)
  // .def("get_molecule", (const RDKit::RWMol &(GemmiSource::*)() const) &
  //                          GemmiSource::get_molecule);

  py::class_<RingData>(_lahuta, "RingData")
      .def(py::init<>())
      .def_readwrite("atom_ids", &RingData::atom_ids)
      .def_property_readonly(
          "center", [](RingData &rd) { return point3d_to_pyarray(rd.center); })
      .def_property_readonly(
          "norm1", [](RingData &rd) { return point3d_to_pyarray(rd.norm1); })
      .def_property_readonly(
          "norm2", [](RingData &rd) { return point3d_to_pyarray(rd.norm2); });

  py::class_<RingDataVec>(_lahuta, "RingDataVec")
      .def(py::init<>())
      .def_readwrite("rings", &RingDataVec::rings)
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
      .def_property_readonly("norm1",
                             [](RingDataVec &rdv) {
                               ssize_t n = rdv.rings.size();
                               ssize_t dim = 3;
                               auto result = py::array_t<double>({n, dim});
                               auto buf = result.request();
                               double *ptr = static_cast<double *>(buf.ptr);

                               for (size_t i = 0; i < n; ++i) {
                                 ptr[i * dim] = rdv.rings[i].norm1.x;
                                 ptr[i * dim + 1] = rdv.rings[i].norm1.y;
                                 ptr[i * dim + 2] = rdv.rings[i].norm1.z;
                               }

                               return result;
                             })
      .def_property_readonly("norm2", [](RingDataVec &rdv) {
        ssize_t n = rdv.rings.size();
        ssize_t dim = 3;
        auto result = py::array_t<double>({n, dim});
        auto buf = result.request();
        double *ptr = static_cast<double *>(buf.ptr);

        for (size_t i = 0; i < n; ++i) {
          ptr[i * dim] = rdv.rings[i].norm2.x;
          ptr[i * dim + 1] = rdv.rings[i].norm2.y;
          ptr[i * dim + 2] = rdv.rings[i].norm2.z;
        }

        return result;
      });

  py::class_<NSResults>(_lahuta, "NSResults")
      .def(py::init<>())
      .def(py::init<const NSResults &>())
      .def(py::init<Pairs &&, std::vector<float> &&>())
      .def(py::init<Pairs &, std::vector<float> &>())
      .def(py::init(
          [](class Luni &luni, Pairs &pairs, std::vector<float> &dists) {
            return NSResults(luni, pairs, dists);
          }))
      .def(py::init(
          [](class Luni &luni, Pairs &&pairs, std::vector<float> &&dists) {
            return NSResults(luni, std::move(pairs), std::move(dists));
          }))
      .def(py::init([](py::list pairs, py::list dists) {
        if (pairs.size() != dists.size()) {
          throw std::invalid_argument(
              "Number of pairs must match number of distances");
        }
        std::vector<std::pair<int, int>> cpp_pairs;
        std::vector<float> cpp_dists;

        for (const auto &pair : pairs) {
          cpp_pairs.emplace_back(pair.cast<std::pair<int, int>>());
        }
        for (const auto &dist : dists) {
          cpp_dists.push_back(dist.cast<float>());
        }

        return NSResults(std::move(cpp_pairs), std::move(cpp_dists));
      }))

      .def("add_neighbors", &NSResults::add_neighbors)
      .def("reserve_space", &NSResults::reserve_space)
      .def("size", &NSResults::size)
      .def("filter", &NSResults::filter)
      // .def("type_filter", &NSResults::type_filter)
      .def("type_filter",
           [](NSResults &nsr, AtomType type, int partner) {
             return nsr.type_filter(type, partner);
           })
      // .def("type_filter2", [](NSResults &nsr, AtomType type, int partner) {
      //   return nsr.type_filter(type, partner);
      // })
      .def("remove_adjascent_pairs",
           &NSResults::remove_adjascent_residueid_pairs)
      .def("clear", &NSResults::clear)
      // .def("get_pairs", &NSResults::get_pairs)
      .def("get_pairs",
           [](NSResults &nsr) {
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
      .def("get_distances_sq", &NSResults::get_distances)
      .def("get_distances",
           [](NSResults &nsr) {
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
      .def("get_luni", &NSResults::get_luni);

  // py::enum_<ContactType>(_lahuta, "ContactType")
  //     .value("AtomAtom", ContactType::AtomAtom)
  //     .value("AtomRing", ContactType::AtomRing);

  py::class_<AtomAtomPair>(_lahuta, "AtomAtomPair")
      .def(py::init<int, int, float>())
      .def_readonly("i", &AtomAtomPair::i)
      .def_readonly("j", &AtomAtomPair::j)
      .def_readonly("d", &AtomAtomPair::d)
      .def("get_pair", &AtomAtomPair::get_pair)
      .def("get_i", &AtomAtomPair::get_i)
      .def("get_j", &AtomAtomPair::get_j)
      .def("names", &AtomAtomPair::names);

  py::class_<AtomRingPair>(_lahuta, "AtomRingPair")
      .def(py::init<int, int, float>())
      .def_readonly("i", &AtomRingPair::i)
      .def_readonly("j", &AtomRingPair::j)
      .def_readonly("d", &AtomRingPair::d)
      .def("get_pair", &AtomRingPair::get_pair)
      .def("get_i", &AtomRingPair::get_i)
      .def("get_j", &AtomRingPair::get_j)
      .def("names", &AtomRingPair::names);

  py::class_<Neighbors<AtomAtomPair>>(_lahuta, "AtomAtomNeighbors")
      .def(py::init([](class Luni &luni, std::vector<AtomAtomPair> data,
                       const AtomAtomPair::RefType &ctx,
                       bool is_sorted = false) {
        return Neighbors<AtomAtomPair>(luni, data, is_sorted);
      }))
      .def(py::init([](class Luni &luni, Pairs pairs, Distances dists,
                       const AtomAtomPair::RefType &ctx,
                       bool is_sorted = false) {
        return Neighbors<AtomAtomPair>(luni, pairs, dists, is_sorted);
      }))
      // .def(py::init([](class Luni &luni, Pairs pairs, Distances dists,
      //                  ContactType c,
      //                  bool is_sorted = false) {
      //   return Neighbors<AtomAtomPair>(luni, pairs, dists, c, is_sorted);
      // }))
      .def_readonly("data", &Neighbors<AtomAtomPair>::data)
      .def("size", &Neighbors<AtomAtomPair>::size)
      .def("type_filter", &Neighbors<AtomAtomPair>::type_filter)
      // .def("type_filter",
      //      [](class Neighbors<AtomAtomPair> &self, int type, int partner) {
      //        return self.type_filter(type, partner);
      //      })

      // .def("remove_adjascent_pairs",
      //      &NSResults::remove_adjascent_residueid_pairs)
      // .def("clear", &NSResults::clear)
      // .def("get_pairs", &NSResults::get_pairs)
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
      // .def("get_distances_sq", &NSResults::get_distances)
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
      .def("intersection_x", [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) {
        return Neighbors<AtomAtomPair>::intersection(self.data, other.data);
      })
      .def("intersection", [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) {
        return self.intersection(other);
      })
      .def("difference", [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) {
        return self.difference(other);
      })
      .def("union", [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) {
        return self.union_(other);
      })
      .def("symmetric_difference", [](Neighbors<AtomAtomPair> &self, Neighbors<AtomAtomPair> &other) {
        return self.symmetric_difference(other);
      })
      .def("get_luni", &Neighbors<AtomAtomPair>::get_luni);


  py::class_<Neighbors<AtomRingPair>>(_lahuta, "AtomRingNeighbors")
      .def(py::init([](class Luni &luni, std::vector<AtomRingPair> data,
                       const AtomRingPair::RefType &ctx,
                       bool is_sorted = false) {
        return Neighbors<AtomRingPair>(luni, data, is_sorted);
      }))
      .def(py::init([](class Luni &luni, Pairs pairs, Distances dists,
                       const AtomRingPair::RefType &ctx,
                       bool is_sorted = false) {
        return Neighbors<AtomRingPair>(luni, pairs, dists, is_sorted);
      }))
      // .def(py::init([](class Luni &luni, Pairs pairs, Distances dists,
      //                  ContactType c,
      //                  bool is_sorted = false) {
      //   return Neighbors<AtomAtomPair>(luni, pairs, dists, c, is_sorted);
      // }))
      .def("size", &Neighbors<AtomRingPair>::size)
      .def("type_filter", &Neighbors<AtomRingPair>::type_filter)
      // .def("type_filter",
      //      [](class Neighbors<AtomAtomPair> &self, int type, int partner) {
      //        return self.type_filter(type, partner);
      //      })

      // .def("remove_adjascent_pairs",
      //      &NSResults::remove_adjascent_residueid_pairs)
      // .def("clear", &NSResults::clear)
      // .def("get_pairs", &NSResults::get_pairs)


  // Neighbors<T> intersection(const Neighbors<T> &other) const {
  //   std::vector<T> result = intersection(data, other.data);
  //   return Neighbors<T>(*m_luni, result);
  // }

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
      // .def("get_distances_sq", &NSResults::get_distances)
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
      .def("get_luni", &Neighbors<AtomRingPair>::get_luni);



  Luni.def(py::init<std::string>())
      .def("find_neighbors_aa", &Luni::find_neighbors<AtomAtomPair>)
      .def("find_neighbors_ar", &Luni::find_neighbors<AtomRingPair>)
      .def("find_neighbors", [](class Luni &luni, float cutoff) {
        return luni._find_neighbors(cutoff);
      }) 
      // .def("find_find_neighbors", &Luni::find_find_neighbors)
      .def("get_atom_types", &Luni::get_atom_types)
      .def("get_rings", &Luni::get_rings)
      .def("match_smarts_string", &Luni::match_smarts_string)
      // .def("atoms", &Luni::atoms,
      // py::return_value_policy::reference_internal)
      .def("n_atoms", &Luni::n_atoms)

      .def(
          "coordinates",
          [](class Luni &luni) {
            return coordinates(
                luni.get_molecule().getConformer().getPositions());
          },
          "Return the coordinates of the molecule as a numpy array")

      .def("names", &Luni::names)
      // .def("names", [](class Luni &luni) {
      //     auto names = luni.names();
      //     py::list result;
      //
      //     for (const auto& name : names) {
      //         result.append(py::str(name));
      //     }
      //
      //     return py::array(result);
      // })
      .def("symbols", &Luni::symbols)
      .def("indices", &Luni::indices)
      .def("elements", &Luni::elements)
      .def("resnames", &Luni::resnames)
      .def("resids", &Luni::resids)
      .def("resindices", &Luni::resindices)
      .def("chainlabels", &Luni::chainlabels)
      .def("get_cutoff", &Luni::get_cutoff);

  py::class_<RDGeom::Point3D>(_lahuta, "Point3D")
      .def(py::init<>())
      .def_readwrite("x", &RDGeom::Point3D::x)
      .def_readwrite("y", &RDGeom::Point3D::y)
      .def_readwrite("z", &RDGeom::Point3D::z);

  // _lahuta.def("atom_type_to_string", &atom_type_to_string);
}

PYBIND11_MODULE(_lahuta, m) {
  m.doc() = "lahuta: A Python binding for the Lahuta library";

  py::class_<RDKit::RWMol> lahutaRWMol(m, "RWMol");
  py::class_<RDKit::Conformer> lahutaConformer(m, "Conformer");

  xbind_atom_types(m);
  bind(m);
}
