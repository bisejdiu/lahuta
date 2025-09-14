#pragma once

#include "compute/parameters.hpp"

// clang-format off
namespace lahuta::topology {
using namespace compute;

namespace param_ids {
  constexpr ParameterInterface::TypeId NEIGHBOR_SEARCH     = 1;
  constexpr ParameterInterface::TypeId RESIDUE_COMPUTATION = 2;
  constexpr ParameterInterface::TypeId BOND_COMPUTATION    = 3;
  constexpr ParameterInterface::TypeId ATOM_TYPING         = 4;
  constexpr ParameterInterface::TypeId RING_COMPUTATION    = 5;
  constexpr ParameterInterface::TypeId NONSTANDARD_BOND_COMPUTATION = 6;
  constexpr ParameterInterface::TypeId SEED_FROM_MODEL     = 7;
}

struct NeighborSearchParams : public ParameterBase<NeighborSearchParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::NEIGHBOR_SEARCH;
  double cutoff = 4.5;
};

struct BondComputationParams : public ParameterBase<BondComputationParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::BOND_COMPUTATION;
};

struct NonStandardBondComputationParams : public ParameterBase<NonStandardBondComputationParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::NONSTANDARD_BOND_COMPUTATION;
};

struct ResidueComputationParams : public ParameterBase<ResidueComputationParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::RESIDUE_COMPUTATION;
};

struct RingComputationParams : public ParameterBase<RingComputationParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::RING_COMPUTATION;
};

struct AtomTypingParams : public ParameterBase<AtomTypingParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::ATOM_TYPING;
  bool use_molstar = true;
};

struct SeedFromModelParams : public ParameterBase<SeedFromModelParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::SEED_FROM_MODEL;
};

} // namespace lahuta::topology
