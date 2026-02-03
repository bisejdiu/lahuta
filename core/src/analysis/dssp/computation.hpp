/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   std::string result;
 *   for (auto part : parts) {
 *     auto* bytes = reinterpret_cast<const std::byte*>(part.data());
 *     for (std::size_t i = 0; i < part.size(); ++i) {
 *       result += static_cast<char>(bytes[i]);
 *     }
 *   }
 *   return result;
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_DSSP_COMPUTATION_HPP
#define LAHUTA_ANALYSIS_DSSP_COMPUTATION_HPP

#include "analysis/dssp/kernel.hpp"
#include "analysis/topology/computation.hpp"
#include "compute/dependency.hpp"
#include "pipeline/task/compute/dynamic_computation.hpp"

namespace lahuta::analysis {
namespace C = lahuta::compute;
namespace P = lahuta::pipeline;

class DsspComputation : public P::DynamicLabelComputation<P::DsspParams, DsspKernel, DsspComputation> {
public:
  using Base = P::DynamicLabelComputation<P::DsspParams, DsspKernel, DsspComputation>;
  using Base::DynamicLabelComputation;

  using dependencies = C::Dependencies<C::Dependency<BuildTopologyComputation, void>>;

  P::DataFieldSet data_requirements() const override {
    return P::DataFieldSet::of({P::DataField::Sequence, P::DataField::Positions});
  }
};

class DsspModelComputation
    : public P::DynamicLabelComputation<P::DsspParams, DsspKernel, DsspModelComputation> {
public:
  using Base = P::DynamicLabelComputation<P::DsspParams, DsspKernel, DsspModelComputation>;
  using Base::DynamicLabelComputation;

  P::DataFieldSet data_requirements() const override {
    return P::DataFieldSet::of({P::DataField::Sequence, P::DataField::Positions});
  }
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_DSSP_COMPUTATION_HPP
