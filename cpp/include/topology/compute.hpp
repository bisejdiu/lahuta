#pragma once

#include "compute/compute_impl.hpp"
#include "compute/dependency.hpp"
#include "nsgrid.hpp"
#include "topology/data.hpp"
#include "topology/kernels.hpp"

// clang-format off
namespace lahuta::topology {
using namespace compute;

template <typename DataT = TopologyData>
class NeighborSearchComputation : public KernelizedROComputation<
    DataT,
    NeighborSearchParams,
    NeighborSearchKernel,
    NeighborSearchComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"neighbors"};
    using dependencies = NoDependencies;
    using NeighborSearchComputation<DataT>::KernelizedROComputation::KernelizedROComputation;
};

template <typename DataT = TopologyData>
class BondComputation : public KernelizedRWComputation<
    DataT,
    BondComputationParams,
    BondKernel,
    BondComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"bonds"};
    using dependencies = Dependencies<Dependency<NeighborSearchComputation<DataT>, std::shared_ptr<lahuta::NSResults>>>;
    using BondComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
};

template <typename DataT = TopologyData>
class ResidueComputation : public KernelizedRWComputation<
    DataT,
    ResidueComputationParams,
    ResidueKernel,
    ResidueComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"residues"};
    using dependencies = NoDependencies;
    using ResidueComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
};

template <typename DataT = TopologyData>
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

template <typename DataT = TopologyData>
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

} // namespace lahuta::topology
