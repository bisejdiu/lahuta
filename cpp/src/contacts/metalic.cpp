#include "contacts/metalic.hpp"
#include "entities/find_contacts.hpp"

namespace lahuta {

bool is_metalic(AtomType at1, AtomType at2) {
  using AtomTypeFlags::has;
  if (has(at1, AtomType::TransitionMetal)) return has(at2, AtomType::DativeBondPartner);
  if (has(at1, AtomType::IonicTypeMetal))  return has(at2, AtomType::IonicTypePartner);
  return false;
}

ContactSet find_metalic(const Topology &topology, const MetalicParams &opts) {
  return find_contacts(
    topology,
    [](const AtomRec &rec) { return (rec.type & (AtomType::IonicTypeMetal   | AtomType::TransitionMetal))   != AtomType::None; },
    [](const AtomRec &rec) { return (rec.type & (AtomType::IonicTypePartner | AtomType::DativeBondPartner)) != AtomType::None; },
    {opts.distance_max, 0, 0, 0.7},
    [&topology, &opts](std::uint32_t idx1, std::uint32_t idx2, float dist) -> InteractionType {
      const auto &m  = topology.atom(idx1);
      const auto &mb = topology.atom(idx2);

      if (dist < opts.distance_max) return InteractionType::None;
      if (idx1 == idx2) return InteractionType::None;

      if (!is_metalic(m.type, mb.type) && !is_metalic(mb.type, m.type)) return InteractionType::None;
      if (topology.molecule().getBondBetweenAtoms(m.idx, mb.idx))       return InteractionType::None;

      return InteractionType::MetalCoordination;
    }
  );
}

} // namespace lahuta
