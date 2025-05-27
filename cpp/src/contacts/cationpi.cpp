#include "contacts/cationpi.hpp"
#include "chemistry/geometry.hpp"
#include "contacts/utils.hpp"
#include "entities/find_contacts.hpp"

// clang-format off
namespace lahuta {

ContactSet find_cationpi(const Topology& topology, const CationPiParams& params) {
  return find_contacts(
    topology,
    [](const GroupRec& rec) { return (rec.a_type & AtomType::POS_IONISABLE) == AtomType::POS_IONISABLE; },
    [](const RingRec & rec)  { return true; }, // we miss positive hits bc we miss genuine aromatic rings in our perception routine
    {params.distance_max, 0.1, 0.5, 1},
    [&topology, &params](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist) -> InteractionType {
      const auto &cation_rec = topology.group(rec_idx_a);
      const auto &ring_rec   = topology.ring(rec_idx_b);

      const auto &mol = topology.molecule();
      const auto *cation_atom = mol.getAtomWithIdx(cation_rec.atoms.front());
      const auto *ring_atom   = mol.getAtomWithIdx(ring_rec.atoms.front());

      if (is_same_residue(mol, *ring_atom, *cation_atom)) return InteractionType::None;

      auto offset = chemistry::compute_in_plane_offset(cation_rec.center, ring_rec.center, ring_rec.normal);
      if (offset > params.offset_max) return InteractionType::None;
      return InteractionType::CationPi;
    }
  );
}


} // namespace lahuta
