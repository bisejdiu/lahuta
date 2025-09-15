#pragma once

#include "compute/compute_impl.hpp"
#include "compute/dependency.hpp"
#include "nsgrid.hpp"
#include "topology/context.hpp"
#include "topology/kernels.hpp"
#include "bonds.hpp"

// clang-format off
namespace lahuta::topology {
using namespace compute;

template <typename DataT = TopologyContext>
class NeighborSearchComputation : public KernelizedROComputation<
    DataT,
    NeighborSearchParams,
    NeighborSearchKernel,
    NeighborSearchComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"neighbors"};
    using dependencies = UnitComputation;
    using NeighborSearchComputation<DataT>::KernelizedROComputation::KernelizedROComputation;
};

template <typename DataT = TopologyContext>
class ResidueComputation : public KernelizedROComputation<
    DataT,
    ResidueComputationParams,
    ResidueKernel,
    ResidueComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"residues"};
    using dependencies = UnitComputation;
    using ResidueComputation<DataT>::KernelizedROComputation::KernelizedROComputation;
};

template <typename DataT = TopologyContext>
class BondComputation : public KernelizedROComputation<
    DataT,
    BondComputationParams,
    BondKernel,
    BondComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"bonds"};
    using dependencies = Dependencies<Dependency<NeighborSearchComputation<DataT>, std::shared_ptr<lahuta::NSResults>>>;
    using BondComputation<DataT>::KernelizedROComputation::KernelizedROComputation;
};

template <typename DataT = TopologyContext>
class NonStandardBondComputation : public KernelizedROComputation<
    DataT,
    NonStandardBondComputationParams,
    NonStandardBondKernel,
    NonStandardBondComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"nonstandard_bonds"};
    using dependencies = Dependencies<Dependency<BondComputation<DataT>, lahuta::BondAssignmentResult>>;
    using NonStandardBondComputation<DataT>::KernelizedROComputation::KernelizedROComputation;
};

template <typename DataT = TopologyContext>
class RingComputation : public KernelizedRWComputation<
    DataT,
    RingComputationParams,
    RingKernel,
    RingComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"rings"};
    using dependencies = Dependencies<Dependency<BondComputation<DataT>,    bool>,
                                      Dependency<ResidueComputation<DataT>, bool>>;
    using RingComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
};

template <typename DataT = TopologyContext>
class AtomTypingComputation : public KernelizedRWComputation<
    DataT,
    AtomTypingParams,
    AtomTypingKernel,
    AtomTypingComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"atom_typing"};
    using dependencies = Dependencies<Dependency<RingComputation<DataT>, bool>>;
    using AtomTypingComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
};

template <typename DataT = TopologyContext>
class SeedFromModelComputation : public KernelizedRWComputation<
    DataT,
    SeedFromModelParams,
    SeedFromModelKernel,
    SeedFromModelComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"seed_from_model"};
    using dependencies = Dependencies<Dependency<ResidueComputation<DataT>, bool>>;
    using SeedFromModelComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
};

template <typename DataT = TopologyContext>
class ModelTopologyComputation : public KernelizedROComputation<
    DataT,
    ModelTopologyParams,
    NoopKernel,
    ModelTopologyComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"model_topology"};
    using dependencies = Dependencies<
        Dependency<ResidueComputation<DataT>, bool>,
        Dependency<SeedFromModelComputation<DataT>, bool>
    >;
    using ModelTopologyComputation<DataT>::KernelizedROComputation::KernelizedROComputation;
};

} // namespace lahuta::topology
