#include "entities.hpp"
#include "atom_types.hpp"
#include "lahuta.hpp"

namespace lahuta {

void AtomEntityCollection::add_data(const RDKit::RWMol &mol, const RDKit::Atom *atom, AtomType type) {
  data.push_back(AtomEntity(type, {atom}, &mol.getConformer().getAtomPos(atom->getIdx()), atom->getIdx()));
}

void RingEntityCollection::add_data(const RDKit::RWMol &mol, const std::vector<int> &ring, int id) {
  RDGeom::Point3D center, norm;
  RingProps::compute_center(&mol, ring, center);
  RingProps::compute_normal(&mol, ring, center, norm);
  std::vector<const RDKit::Atom *> atoms;
  for (const int &atom_idx : ring) {
    atoms.push_back(mol.getAtomWithIdx(atom_idx));
  }
  data.push_back(RingEntity(center, norm, atoms, id));
}

void GroupEntityCollection::add_data(
    AtomType type, FeatureGroup group, const std::vector<const RDKit::Atom *> &members) {
  auto center = RDGeom::Point3D(0.0, 0.0, 0.0);
  for (const auto *atom : members) {
    center += atom->getOwningMol().getConformer().getAtomPos(atom->getIdx());
  }
  data.emplace_back(type, group, members, center / members.size());
}

// overwrite the default implementation, because we store pointers to atom positions
const RDGeom::POINT3D_VECT AtomEntityCollection::positions() const {
  RDGeom::POINT3D_VECT pos;
  pos.reserve(data.size());
  for (const auto &atom_data : data) {
    pos.push_back(*atom_data.center);
  }
  return pos;
}

/*std::vector<std::vector<double>> RingDataVec::norm() const {*/
/*  std::vector<std::vector<double>> normals;*/
/*  for (const auto &ring : rings) {*/
/*    normals.push_back({ring.norm.x, ring.norm.y, ring.norm.z});*/
/*  }*/
/*  return normals;*/
/*}*/

std::vector<double> RingEntityCollection::compute_angles(
    const std::vector<int> &ring_indices, const std::vector<std::vector<double>> &points) const {
  std::vector<double> angles;
  angles.reserve(ring_indices.size());
  for (size_t i = 0; i < ring_indices.size(); ++i) {
    auto angle = compute_angle(data[ring_indices[i]], points[i]);
    angles.push_back(angle);
  }
  return angles;
}

double RingEntityCollection::compute_angle(const RingEntity &rd, const std::vector<double> &_point) const {
  RDGeom::Point3D point(_point[0], _point[1], _point[2]);
  auto vector_point_to_plane = point - rd.center;
  vector_point_to_plane.normalize();

  double cos_theta = vector_point_to_plane.dotProduct(rd.norm);
  cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

  double theta_radians = std::acos(cos_theta); // in radians
  return theta_radians * (180.0 / M_PI);
}

const AtomEntityCollection AtomEntityCollection::filter(const Luni *luni, AtomType type, FeatureTypeCheckFunc check_func) {
  const std::vector<AtomType> &atom_types = luni->get_atom_types();
  const auto &mol = luni->get_molecule();

  AtomEntityCollection atom_data_vec;
  for (const auto *atom : mol.atoms()) {
    if (check_func(atom_types[atom->getIdx()], type)) {
      atom_data_vec.add_data(mol, atom, atom_types[atom->getIdx()]);
    }
  }
  return std::move(atom_data_vec);
}

const GroupEntityCollection GroupEntityCollection::filter(const Luni* luni, AtomType type, FeatureTypeCheckFunc check_func) {
  const GroupEntityCollection &features = luni->get_features();
  GroupEntityCollection feature_vec;
  for (const auto &feature : features.get_data()) {
    if (check_func(feature.type, type)) {
      feature_vec.add_data(feature);
    }
  }

  return feature_vec;
}


RingEntityCollection get_rings(const Luni *luni) { return luni->get_rings(); }

std::vector<const RDKit::Atom *> get_atom_types(const Luni *luni, AtomType type) {
  const std::vector<AtomType> &atom_types = luni->get_atom_types();
  const auto &mol = luni->get_molecule();

  std::vector<const RDKit::Atom *> filtered_atoms;
  for (const auto *atom : mol.atoms()) {
    if (AtomTypeFlags::has_any(atom_types[atom->getIdx()], type)) {
      filtered_atoms.push_back(atom);
    }
  }
  return filtered_atoms;
}

} // namespace lahuta
