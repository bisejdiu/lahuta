#include "contacts/pistacking.hpp"
#include "chemistry/geometry.hpp"
#include "contacts/utils.hpp"
#include "entities/find_contacts.hpp"

// clang-format off
namespace lahuta {

ContactSet find_pistacking(const Topology &topology, const PiStackingParams& params) {
  return find_contacts(
    topology,
    [](const RingRec& rec) { return rec.aromatic; },
    {params.distance_max, 0.1, 1},
    [&topology, &params](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist) -> InteractionType {

      const auto& ring_rec_a = topology.ring(rec_idx_a);
      const auto& ring_rec_b = topology.ring(rec_idx_b);

      const auto& mol = topology.molecule();
      const auto* r_atom_a = mol.getAtomWithIdx(ring_rec_a.atoms.front());
      const auto* r_atom_b = mol.getAtomWithIdx(ring_rec_b.atoms.front());

      // e.g. trp
      if (is_same_residue(mol, *r_atom_a, *r_atom_b)) return InteractionType::None;

      auto dot_product = ring_rec_a.normal.dotProduct(ring_rec_b.normal);
      auto angle = std::acos(std::clamp(dot_product, -1.0, 1.0));
      if (angle > M_PI / 2) angle = M_PI - angle; // obtuse -> acute

      double offset_a = chemistry::compute_in_plane_offset(ring_rec_a.center, ring_rec_b.center, ring_rec_a.normal);
      double offset_b = chemistry::compute_in_plane_offset(ring_rec_b.center, ring_rec_a.center, ring_rec_b.normal);

      if (std::min(offset_a, offset_b) > params.offset_max) return InteractionType::None;
      if (angle <= params.angle_dev_max)                    return InteractionType::PiStackingP;
      if (std::abs(angle - M_PI/2) <= params.angle_dev_max) return InteractionType::PiStackingT;

      return InteractionType::None;
    }
  );
}


} // namespace lahuta
