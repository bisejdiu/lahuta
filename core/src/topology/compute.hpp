/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto parts = std::make_tuple("besian", "sejdiu", "@gmail.com");
 *   return std::apply([](auto... p) { std::string s; (s.append(p), ...); return s; }, parts);
 * }();
 *
 */

#ifndef LAHUTA_TOPOLOGY_COMPUTE_HPP
#define LAHUTA_TOPOLOGY_COMPUTE_HPP

#include "bonds/bonds.hpp"
#include "compute/compute_impl.hpp"
#include "compute/dependency.hpp"
#include "spatial/nsresults.hpp"
#include "topology/context.hpp"
#include "topology/kernels.hpp"

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
class ResidueComputation : public KernelizedRWComputation<
    DataT,
    ResidueComputationParams,
    ResidueKernel,
    ResidueComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"residues"};
    using dependencies = UnitComputation;
    using ResidueComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
};

template <typename DataT = TopologyContext>
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

template <typename DataT = TopologyContext>
class NonStandardBondComputation : public KernelizedRWComputation<
    DataT,
    NonStandardBondComputationParams,
    NonStandardBondKernel,
    NonStandardBondComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"nonstandard_bonds"};
    using dependencies = Dependencies<Dependency<BondComputation<DataT>, lahuta::BondAssignmentResult>>;
    using NonStandardBondComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
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

#endif // LAHUTA_TOPOLOGY_COMPUTE_HPP
