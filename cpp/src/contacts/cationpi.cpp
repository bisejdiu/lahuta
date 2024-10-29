#include "lahuta.hpp"
#include "nn.hpp"
#include "contacts/cationpi.hpp"

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

void find_cationpi(const Luni *luni, GeometryOptions opts, Contacts &contacts) {
  const auto &mol = luni->get_molecule();
  const auto rings = luni->get_rings();

  auto features = get_features(luni, AtomType::POS_IONISABLE);
  if (features.features.empty()) {
    return;
  }

  double cationpi_max_dist = 6.0;
  auto grid = FastNS(rings.centers(), cationpi_max_dist);
  auto nbrs = grid.search(features.positions());

  for (const auto &[pair, dist] : nbrs) {
    auto [feature_index, ring_index] = pair;
    const auto &ring = rings.rings[ring_index];
    const auto &feature = features[feature_index];

    auto first_ring_atom = ring.atoms.front();

    // FIX: does this make sense here?
    if (is_same_residue(mol, *first_ring_atom, *feature.members.front())) {
      continue;
    }

    auto offset = compute_in_plane_offset(feature.center, ring.center, ring.norm);

    if (offset <= 2.2) {
      EntityID entity1 = make_entity_id(EntityType::Group, feature.get_id());
      EntityID entity2 = make_entity_id(EntityType::Ring, ring.get_id());

      contacts.add(Contact(
        entity1,
        entity2,
        dist,
        InteractionType::Ionic));
    }
  }
}

} // namespace lahuta
