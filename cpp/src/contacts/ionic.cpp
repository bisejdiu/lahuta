#include "contacts/ionic.hpp"
#include "entities/find_contacts.hpp"

// clang-format off
namespace lahuta {

ContactSet find_ionic(const Topology& topology, const IonicParams& params) {
  return find_contacts(
    topology,
    [](const GroupRec &r) { return (r.a_type & AtomType::POS_IONISABLE) == AtomType::POS_IONISABLE; },
    [](const GroupRec &r) { return (r.a_type & AtomType::NEG_IONISABLE) == AtomType::NEG_IONISABLE; },
    {params.distance_max, 0.5, 0.5},
    [&topology, &params](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist) -> InteractionType {
      if (dist < 2.0 || rec_idx_a == rec_idx_b) return InteractionType::None; // FIX: idx1 == idx2 should not be necessary
      return InteractionType::Ionic;
    }
  );
}

} // namespace lahuta
