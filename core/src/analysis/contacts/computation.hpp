#pragma once

#include "analysis/contacts/kernel.hpp"
#include "analysis/topology/computation.hpp"
#include "compute/dependency.hpp"
#include "pipeline/compute/dynamic_computation.hpp"
#include "pipeline/compute/parameters.hpp"

// clang-format off
namespace lahuta::analysis::contacts {
using namespace lahuta::topology::compute;

// Uses dynamic labels but has static dependency on topology computation.
class ContactsComputation : public DynamicLabelComputation<ContactsParams, ContactsKernel, ContactsComputation> {
public:
  using Base = DynamicLabelComputation<ContactsParams, ContactsKernel, ContactsComputation>;
  using Base::DynamicLabelComputation;
  using dependencies = Dependencies<Dependency<analysis::topology::BuildTopologyComputation, void>>;
};

} // namespace lahuta::analysis::contacts
