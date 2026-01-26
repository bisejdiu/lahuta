#ifndef LAHUTA_ANALYSIS_CONTACTS_COMPUTATION_HPP
#define LAHUTA_ANALYSIS_CONTACTS_COMPUTATION_HPP

#include "analysis/contacts/kernel.hpp"
#include "analysis/topology/computation.hpp"
#include "compute/dependency.hpp"
#include "pipeline/task/compute/dynamic_computation.hpp"
#include "pipeline/task/compute/parameters.hpp"

namespace lahuta::analysis {
namespace C = lahuta::compute;
namespace P = lahuta::pipeline;

// Uses dynamic labels but has static dependency on topology computation.
// clang-format off
class ContactsComputation : public P::DynamicLabelComputation<P::ContactsParams, ContactsKernel, ContactsComputation> {
public:
  using Base = P::DynamicLabelComputation<P::ContactsParams, ContactsKernel, ContactsComputation>;
  using Base::DynamicLabelComputation;
  using dependencies = C::Dependencies<C::Dependency<BuildTopologyComputation, void>>;
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_CONTACTS_COMPUTATION_HPP
