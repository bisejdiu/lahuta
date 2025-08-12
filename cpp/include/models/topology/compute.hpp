#ifndef LAHUTA_MODELS_TOPOLOGY_COMPUTE_HPP
#define LAHUTA_MODELS_TOPOLOGY_COMPUTE_HPP

#include "compute/compute_impl.hpp"
#include "compute/dependency.hpp"
#include "models/topology/data.hpp"
#include "models/topology/kernels.hpp"
#include "models/topology/parameters.hpp"

// clang-format off
namespace lahuta::models::topology {

template <typename DataT = ModelData>
class ModelAtomsComputation : public KernelizedRWComputation<
    DataT,
    ModelAtomsParams,
    ModelAtomsKernel,
    ModelAtomsComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"model_atoms"};
    using dependencies = UnitComputation;
    using ModelAtomsComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
};

template <typename DataT = ModelData>
class ModelBondsComputation : public KernelizedRWComputation<
    DataT,
    ModelBondsParams,
    ModelBondsKernel,
    ModelBondsComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"model_bonds"};
    using dependencies = Dependencies<Dependency<ModelAtomsComputation<DataT>, bool>>;
    using ModelBondsComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
};

template <typename DataT = ModelData>
class ModelPositionsComputation : public KernelizedRWComputation<
    DataT,
    ModelPositionsParams,
    ModelPositionsKernel,
    ModelPositionsComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"model_positions"};
    using dependencies = Dependencies<Dependency<ModelAtomsComputation<DataT>, bool>>;
    using ModelPositionsComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
};

template <typename DataT = ModelData>
class ModelAromaticsComputation : public KernelizedRWComputation<
    DataT,
    ModelAromaticsParams,
    ModelAromaticsKernel,
    ModelAromaticsComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"model_aromatics"};
    using dependencies = Dependencies<Dependency<ModelBondsComputation<DataT>, bool>>;
    using ModelAromaticsComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
};

template <typename DataT = ModelData>
class ModelDisulfidesComputation : public KernelizedRWComputation<
    DataT,
    ModelDisulfidesParams,
    ModelDisulfidesKernel,
    ModelDisulfidesComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"model_disulfides"};
    using dependencies = Dependencies<Dependency<ModelBondsComputation<DataT>, bool>,
                                      Dependency<ModelPositionsComputation<DataT>, bool>>;
    using ModelDisulfidesComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
};

template <typename DataT = ModelData>
class ModelBuildComputation : public KernelizedRWComputation<
    DataT,
    ModelBuildParams,
    ModelBuildKernel,
    ModelBuildComputation<DataT>> {
public:
    constexpr static const ComputationLabel label{"model_build"};
    using dependencies = Dependencies<Dependency<ModelAtomsComputation<DataT>, bool>>;
    using ModelBuildComputation<DataT>::KernelizedRWComputation::KernelizedRWComputation;
};

} // namespace lahuta::models::topology

#endif // LAHUTA_MODELS_TOPOLOGY_COMPUTE_HPP
