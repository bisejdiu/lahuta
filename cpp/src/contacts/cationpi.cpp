#include "contacts/cationpi.hpp"
#include "contacts/search.hpp"
#include "lahuta.hpp"
#include "nn.hpp"

namespace lahuta {

bool is_same_residue(const RDKit::RWMol &mol, const RDKit::Atom &atom_a, const RDKit::Atom &atom_b) {
  auto info_a = static_cast<const RDKit::AtomPDBResidueInfo *>(atom_a.getMonomerInfo());
  auto info_b = static_cast<const RDKit::AtomPDBResidueInfo *>(atom_b.getMonomerInfo());
  return info_a->getResidueNumber() == info_b->getResidueNumber();
}

// Project a vector onto a plane defined by a normal vector
RDGeom::Point3D project_on_plane(const RDGeom::Point3D &vector, const RDGeom::Point3D &plane_normal) {
  // subtract component along the normal
  double scalar_projection = vector.dotProduct(plane_normal);
  RDGeom::Point3D projected_vec = vector - (plane_normal * scalar_projection);
  return projected_vec;
}

double compute_in_plane_offset(
    const RDGeom::Point3D &pos_a, const RDGeom::Point3D &pos_b, const RDGeom::Point3D &normal) {
  RDGeom::Point3D vec_ab = pos_a - pos_b;
  RDGeom::Point3D projected_vec = project_on_plane(vec_ab, normal);
  double in_plane_offset = projected_vec.length();
  return in_plane_offset;
}

Contacts find_cationpi(const Luni &luni, CationPiParams opts) {

  Contacts contacts(&luni);
  const auto &mol = luni.get_molecule();
  const auto rings = luni.get_rings();

  auto features = get_features(&luni, AtomType::POS_IONISABLE);

  // FIX: early return does not provide important performance gain.
  // We should just let it run, if we can be sure no errors will occur.
  if (features.get_data().empty()) {
    return contacts;
  }

  EntityNeighborSearch ens(mol.getConformer());
  auto nbrs = ens.search(features, rings, opts.distance_max);

  for (const auto &[pair, dist] : nbrs) {
    auto [feature_index, ring_index] = pair;
    const auto &ring = rings.rings[ring_index];
    const auto &feature = features[feature_index];

    auto first_ring_atom = ring.atoms.front();

    // FIX: does this make sense here?
    if (is_same_residue(mol, *first_ring_atom, *feature.members.front())) continue;

    auto offset = compute_in_plane_offset(feature.center, ring.center, ring.norm);

    // FIX: Are we sure we won't get multiple contacts for the same feature?
    if (offset <= opts.offset_max) {
      contacts.add(Contact(
          make_entity_id(EntityType::Group, feature.get_id()),
          make_entity_id(EntityType::Ring, ring.get_id()),
          dist,
          InteractionType::CationPi));
    }
  }

  return contacts;
}

} // namespace lahuta
