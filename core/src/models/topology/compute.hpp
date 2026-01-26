#ifndef LAHUTA_MODELS_TOPOLOGY_COMPUTE_HPP
#define LAHUTA_MODELS_TOPOLOGY_COMPUTE_HPP

#include "compute/compute_impl.hpp"
#include "compute/dependency.hpp"
#include "models/topology/data.hpp"
#include "models/topology/kernels.hpp"
#include "models/topology/parameters.hpp"

namespace lahuta::models::topology {
namespace C = lahuta::compute;

// clang-format off
template <typename DataT = ModelData>
class ModelAtomsComputation : public C::KernelizedRWComputation<
    DataT,
    ModelAtomsParams,
    ModelAtomsKernel,
    ModelAtomsComputation<DataT>> {
public:
    constexpr static const C::ComputationLabel label{"model_atoms"};
    using dependencies = C::UnitComputation;
    using ModelAtomsComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
};

template <typename DataT = ModelData>
class ModelBondsComputation : public C::KernelizedRWComputation<
    DataT,
    ModelBondsParams,
    ModelBondsKernel,
    ModelBondsComputation<DataT>> {
public:
    constexpr static const C::ComputationLabel label{"model_bonds"};
    using dependencies = C::Dependencies<C::Dependency<ModelAtomsComputation<DataT>, bool>>;
    using ModelBondsComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
};

template <typename DataT = ModelData>
class ModelPositionsComputation : public C::KernelizedRWComputation<
    DataT,
    ModelPositionsParams,
    ModelPositionsKernel,
    ModelPositionsComputation<DataT>> {
public:
    constexpr static const C::ComputationLabel label{"model_positions"};
    using dependencies = C::Dependencies<C::Dependency<ModelAtomsComputation<DataT>, bool>>;
    using ModelPositionsComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
};

template <typename DataT = ModelData>
class ModelAromaticsComputation : public C::KernelizedRWComputation<
    DataT,
    ModelAromaticsParams,
    ModelAromaticsKernel,
    ModelAromaticsComputation<DataT>> {
public:
    constexpr static const C::ComputationLabel label{"model_aromatics"};
    using dependencies = C::Dependencies<C::Dependency<ModelBondsComputation<DataT>, bool>>;
    using ModelAromaticsComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
};

template <typename DataT = ModelData>
class ModelDisulfidesComputation : public C::KernelizedRWComputation<
    DataT,
    ModelDisulfidesParams,
    ModelDisulfidesKernel,
    ModelDisulfidesComputation<DataT>> {
public:
    constexpr static const C::ComputationLabel label{"model_disulfides"};
    using dependencies = C::Dependencies<C::Dependency<ModelBondsComputation<DataT>, bool>,
                                      C::Dependency<ModelPositionsComputation<DataT>, bool>>;
    using ModelDisulfidesComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
};

template <typename DataT = ModelData>
class ModelBuildComputation : public C::KernelizedRWComputation<
    DataT,
    ModelBuildParams,
    ModelBuildKernel,
    ModelBuildComputation<DataT>> {
public:
    constexpr static const C::ComputationLabel label{"model_build"};
    using dependencies = C::Dependencies<C::Dependency<ModelAtomsComputation<DataT>, bool>>;
    using ModelBuildComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
};

} // namespace lahuta::models::topology

#endif // LAHUTA_MODELS_TOPOLOGY_COMPUTE_HPP
