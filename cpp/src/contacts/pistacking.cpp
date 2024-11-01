#include "contacts/pistacking.hpp"
#include "contacts/cationpi.hpp"
#include "contacts/search.hpp"
#include "lahuta.hpp"

namespace lahuta {

bool is_same_residue(const RDKit::RWMol &mol, const RingData &ring_a, const RingData &ring_b) {
  auto atom_a = ring_a.atoms[0];
  auto atom_b = ring_b.atoms[0];
  auto info_a = static_cast<const RDKit::AtomPDBResidueInfo *>(atom_a->getMonomerInfo());
  auto info_b = static_cast<const RDKit::AtomPDBResidueInfo *>(atom_b->getMonomerInfo());
  return info_a->getResidueNumber() == info_b->getResidueNumber();
}

void find_pistacking(const Luni *luni, GeometryOptions opts, Contacts &contacts) {

  const auto &mol = luni->get_molecule();
  const auto rings = luni->get_rings();

  double dist_max = 6.0;
  EntityNeighborSearch ens(mol.getConformer());
  auto nbrs = ens.search(rings, dist_max);

  for (const auto &[pair, dist] : nbrs) {
    auto [ring_index_a, ring_index_b] = pair;
    const auto &ring_a = rings[ring_index_a];
    const auto &ring_b = rings[ring_index_b];

    if (is_same_residue(mol, ring_a, ring_b)) { // e.g. TRP
      continue;
    }

    auto dot_product = ring_a.norm.dotProduct(ring_b.norm);
    auto angle = std::acos(std::clamp(dot_product, -1.0, 1.0));
    if (angle > M_PI / 2) angle = M_PI - angle;

    double offset_a = compute_in_plane_offset(ring_a.center, ring_b.center, ring_a.norm);
    double offset_b = compute_in_plane_offset(ring_b.center, ring_a.center, ring_b.norm);

    double offset = std::min(offset_a, offset_b);

    if (offset <= offset_max) {
      EntityID entity1 = make_entity_id(EntityType::Ring, ring_a.get_id());
      EntityID entity2 = make_entity_id(EntityType::Ring, ring_b.get_id());
      if (angle <= AngleDevMax) {
        std::cout << "Found PiStacking: Parallel" << std::endl;
        contacts.add(Contact(entity1, entity2, dist, InteractionType::PiStackingP));
      } else if (std::abs(angle - M_PI / 2) <= AngleDevMax) {
        std::cout << "Found PiStacking: T-Shaped" << std::endl;
        contacts.add(Contact(entity1, entity2, dist, InteractionType::PiStackingT));
      }
    }
  }
}
} // namespace lahuta
