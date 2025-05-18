#ifndef LAHUTA_ENTITIES_INTERACTION_TYPES_HPP
#define LAHUTA_ENTITIES_INTERACTION_TYPES_HPP

namespace lahuta {

enum class InteractionType {
  None,
  Generic,
  Hydrophobic,
  Halogen,
  HydrogenBond,
  WeakHydrogenBond,
  Ionic,
  MetalCoordination,
  CationPi,
  PiStackingP,
  PiStackingT,
};

} // namespace lahuta

#endif // LAHUTA_ENTITIES_INTERACTION_TYPES_HPP
