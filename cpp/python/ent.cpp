#include <pybind11/complex.h> // FIX: is this needed?
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/RDKitBase.h>

#include "typing/types.hpp"
#include "array.hpp"
#include "contacts.hpp"
#include "contacts/search.hpp"
#include "entities.hpp"
#include "lahuta.hpp"

#define BIND_SEARCH_ONE(T)     .def("search", py::overload_cast<const T &,              double>(&EntityNeighborSearch::template search<T>))
#define BIND_SEARCH_TWO(T1, T2).def("search", py::overload_cast<const T1 &, const T2 &, double>(&EntityNeighborSearch::template search<T1, T2>))

using namespace lahuta;
namespace py = pybind11;

template <typename Collection, typename Func>
py::array_t<double> extract_vectors(const Collection &coll, Func extractor) {
  const ssize_t n = coll.data.size();
  const ssize_t dim = 3;
  auto result = py::array_t<double>({n, dim});
  auto buf = result.request();
  double *ptr = static_cast<double *>(buf.ptr);
  for (ssize_t i = 0; i < n; ++i) {
    // extractor returns an object with members x, y, and z.
    auto vec = extractor(coll.data[i]);
    ptr[i * dim] = vec.x;
    ptr[i * dim + 1] = vec.y;
    ptr[i * dim + 2] = vec.z;
  }
  return result;
}

// clang-format off
void bind_entities(py::module &_lahuta) {
  py::enum_<FeatureGroup> FeatureGroup_(_lahuta, "FeatureGroup");

  py::class_<AtomEntity>  AtomEntity_ (_lahuta, "AtomEntity");
  py::class_<RingEntity>  RingEntity_ (_lahuta, "RingData");
  py::class_<GroupEntity> GroupEntity_(_lahuta, "GroupEntity");
  py::class_<AtomEntityCollection>  AtomEntityCollection_ (_lahuta, "AtomEntityCollection");
  py::class_<RingEntityCollection>  RingEntityCollection_ (_lahuta, "RingEntityCollection");
  py::class_<GroupEntityCollection> GroupEntityCollection_(_lahuta, "GroupEntityCollection");
  py::class_<EntityNeighborSearch>  EntityNeighborSearch_ (_lahuta, "EntityNeighborSearch");

  EntityNeighborSearch_
    BIND_SEARCH_TWO(AtomEntityCollection,  AtomEntityCollection)
    BIND_SEARCH_TWO(AtomEntityCollection,  RingEntityCollection)
    BIND_SEARCH_TWO(AtomEntityCollection,  GroupEntityCollection)
    BIND_SEARCH_TWO(RingEntityCollection,  AtomEntityCollection)
    BIND_SEARCH_TWO(RingEntityCollection,  RingEntityCollection)
    BIND_SEARCH_TWO(RingEntityCollection,  GroupEntityCollection)
    BIND_SEARCH_TWO(GroupEntityCollection, AtomEntityCollection)
    BIND_SEARCH_TWO(GroupEntityCollection, RingEntityCollection)
    BIND_SEARCH_TWO(GroupEntityCollection, GroupEntityCollection)

    BIND_SEARCH_ONE(AtomEntityCollection)
    BIND_SEARCH_ONE(RingEntityCollection)
    BIND_SEARCH_ONE(GroupEntityCollection)
  ;

  FeatureGroup_
    .value("None",            FeatureGroup::None)
    .value("QuaternaryAmine", FeatureGroup::QuaternaryAmine)
    .value("TertiaryAmine",   FeatureGroup::TertiaryAmine)
    .value("Sulfonium",       FeatureGroup::Sulfonium)
    .value("SulfonicAcid",    FeatureGroup::SulfonicAcid)
    .value("Sulfate",         FeatureGroup::Sulfate)
    .value("Phosphate",       FeatureGroup::Phosphate)
    .value("Halocarbon",      FeatureGroup::Halocarbon)
    .value("Guanidine",       FeatureGroup::Guanidine)
    .value("Acetamidine",     FeatureGroup::Acetamidine)
    .value("Carboxylate",     FeatureGroup::Carboxylate);


  AtomEntity_
    .def(py::init([](const RDKit::Atom * atom, AtomType at, RDGeom::Point3D *center, size_t id) { return new AtomEntity(at, {atom}, center, id); }))
    .def(py::init([](const class Luni& luni, int atom_idx, AtomType at, RDGeom::Point3D* center, size_t id) { 
      return new AtomEntity(at, {luni.get_atom(atom_idx)}, center, id); 
    }))

    .def_property_readonly("id",     &AtomEntity::get_id)
    .def_property_readonly("center", [](AtomEntity &ae) {return point3d_to_pyarray(*ae.center);})

    .def_readwrite("type", &AtomEntity::type)
    .def_readwrite("atom", &AtomEntity::atom)

    .def("get_data", &AtomEntity::get_data, py::return_value_policy::reference)
    .def("has_atom", &AtomEntity::has_atom)

    .def("__str__", [](const AtomEntity &ae) { return "AtomEntity(ID: " + std::to_string(ae.get_id()) + ", " + atom_type_to_string(ae.type) + ")"; });


  RingEntity_
    .def(py::init([](const std::vector<const RDKit::Atom*> &atoms, RDGeom::Point3D &center, RDGeom::Point3D &normal, size_t id) {
      return new RingEntity(center, normal, atoms, id);
    }))
    .def(py::init([](const class Luni& luni, const std::vector<int> &atom_ids, RDGeom::Point3D &center, RDGeom::Point3D &norm, size_t id) {
      std::vector<const RDKit::Atom *> atoms;
      for (const int &atom_idx : atom_ids) {
        atoms.push_back(luni.get_atom(atom_idx));
      }
      return new RingEntity(center, norm, atoms, id);
    }))

    .def_property_readonly("id",     &RingEntity::get_id)
    .def_property_readonly("center", [](RingEntity &re) {return point3d_to_pyarray(re.center);})
    .def_property_readonly("normal", [](RingEntity &re) {return point3d_to_pyarray(re.normal);})

    /*.def_readwrite("type", &RingEntity::type)*/
    .def_readwrite("atoms", &RingEntity::atoms)

    .def("get_data", &RingEntity::get_data, py::return_value_policy::reference)
    .def("has_atom", &RingEntity::has_atom)

    .def("get_atom_ids", [](RingEntity &rd) {
      std::vector<int> ids;
      ids.reserve(rd.atoms.size());
      for (const auto *atom : rd.atoms) {
        ids.push_back(atom->getIdx());
      }
      return ids;
    })
    .def("__str__", [](const RingEntity &re) { return "RingEntity(" + std::to_string(re.get_id()) + " >ring type not yet supported<)"; });


  GroupEntity_
    .def(py::init([](std::vector<const RDKit::Atom*> &atoms, RDGeom::Point3D &center, AtomType type, FeatureGroup group, size_t id) {
      auto ge = new GroupEntity(type, group, atoms, center);
      ge->set_id(id);
      return ge;
    }))
    .def(py::init([](const class Luni &luni, const std::vector<int> &atom_ids, AtomType type, FeatureGroup group, const RDGeom::Point3D &center, size_t id) {
      std::vector<const RDKit::Atom *> atoms;
      for (const int &atom_idx : atom_ids) {
        atoms.push_back(luni.get_molecule().getAtomWithIdx(atom_idx));
      }
      auto ge = new GroupEntity(type, group, atoms, center);
      ge->set_id(id);
      return ge;
    }))

    .def_property_readonly("id", &GroupEntity::get_id)
    .def_property_readonly("center", [](GroupEntity &ge) {return point3d_to_pyarray(ge.center);})

    .def_readwrite("type",  &GroupEntity::type)
    .def_readwrite("group", &GroupEntity::group)
    .def_readwrite("atoms", &GroupEntity::atoms)

    .def("get_data", &GroupEntity::get_data, py::return_value_policy::reference)
    .def("has_atom", &GroupEntity::has_atom)

    .def("get_atom_ids", [](GroupEntity &ge) {
      std::vector<int> ids;
      ids.reserve(ge.atoms.size());
      for (const auto *atom : ge.atoms) {
        ids.push_back(atom->getIdx());
      }
      return ids;
    })
    .def("__str__", [](const GroupEntity &ge) { return "GroupEntity(ID: " + std::to_string(ge.get_id()) + ", " + atom_type_to_string(ge.type) + ")"; });


  AtomEntityCollection_
    .def(py::init<>())

    .def("add_entity",    [](AtomEntityCollection &aec, AtomEntity &et) { aec.add_data(et);} )
    .def("get_entities",  &AtomEntityCollection::get_data)
    .def("get_size",      &AtomEntityCollection::size)
    .def("get_positions", &AtomEntityCollection::positions) // FIX: should return a NumPy array
    .def("get_atom_ids",  &AtomEntityCollection::atom_ids)
    .def("get_atoms",     &AtomEntityCollection::atoms, py::return_value_policy::reference)

    .def_static("filter", &AtomEntityCollection::filter, py::arg("luni"), py::arg("type"), py::arg("check_func") = static_cast<FeatureTypeCheckFunc>(AtomTypeFlags::has_any));


  RingEntityCollection_
    .def(py::init<>())
    .def("add_entity",    [](RingEntityCollection &rec, RingEntity &re) { rec.add_data(re); })
    .def("get_entities",  &RingEntityCollection::get_data)
    .def("get_size",      &RingEntityCollection::size)
    .def("get_positions", &RingEntityCollection::positions)
    .def("get_atom_ids",  &RingEntityCollection::atom_ids)
    .def("get_atoms",     &RingEntityCollection::atoms, py::return_value_policy::reference)

    .def("compute_angles", &RingEntityCollection::compute_angles, py::arg("ring_indices"), py::arg("points"))
    .def("compute_angle",  &RingEntityCollection::compute_angle,  py::arg("ring_entity"),  py::arg("point"))

    // FIX: Should RingEntityCollection implement `filter`? Could we filter: RingEntityCollection.filter(luni, AtomType.AROMATIC)
    .def_property_readonly("rings",   [](RingEntityCollection &rdc) { return rdc.get_data(); })
    .def_property_readonly("centers", [](RingEntityCollection &rec) { return extract_vectors(rec, [](const auto &e) -> const auto& { return e.center; });})
    .def_property_readonly("normal",  [](RingEntityCollection &rec) { return extract_vectors(rec, [](const auto &e) -> const auto& { return e.normal; });});


  GroupEntityCollection_
    .def(py::init<>())
    .def(py::init<std::vector<GroupEntity>>())

    .def_static("filter", &GroupEntityCollection::filter, py::arg("luni"), py::arg("type"), py::arg("check_func") = static_cast<FeatureTypeCheckFunc>(AtomTypeFlags::has_any))

    .def("add_entity",    [](GroupEntityCollection &gec, GroupEntity &ge) { gec.add_data(ge); })
    .def("get_entities",  &GroupEntityCollection::get_data)
    .def("get_size",      &GroupEntityCollection::size)
    .def("get_positions", &GroupEntityCollection::positions)
    .def("get_atom_ids",  &GroupEntityCollection::atom_ids)
    .def("get_atoms",     &GroupEntityCollection::atoms, py::return_value_policy::reference)

    .def_static("filter", &GroupEntityCollection::filter, py::arg("luni"), py::arg("type"), py::arg("check_func") = static_cast<FeatureTypeCheckFunc>(AtomTypeFlags::has_any))
    .def_property_readonly("centers", [](GroupEntityCollection &gec) { return extract_vectors(gec, [](const auto &entity) -> const auto& { return entity.center; }); });

}

