#pragma once

#include "analysis/contacts/hooks.hpp"
#include "analysis/contacts/kernel.hpp"
#include "analysis/topology/computation.hpp"
#include "compute/dependency.hpp"
#include "pipeline/compute/dynamic_computation.hpp"
#include "pipeline/compute/parameters.hpp"
#include "pipeline/dynamic/types.hpp"

// clang-format off
namespace lahuta::analysis::contacts {
using namespace lahuta::topology::compute;

// Uses dynamic labels but has static dependency on topology computation.
class ContactsComputation : public DynamicLabelComputation<ContactsParams, ContactsKernel, ContactsComputation> {
public:
  using Base = DynamicLabelComputation<ContactsParams, ContactsKernel, ContactsComputation>;
  using Base::DynamicLabelComputation;
  using dependencies = Dependencies<Dependency<analysis::topology::BuildTopologyComputation, void>>;

  // Post-completion hook tied to this analysis.
  void on_complete(DataContext<PipelineContext, Mut::ReadWrite>& ctx, ComputationResult& res) override {
    auto* tctx = ctx.data().ctx;
    if (!tctx) return;
    // Build summary JSON from context keys set by ContactsKernel
    auto json = build_contacts_summary_json(*tctx);
    pipeline::dynamic::Emission e{"summary", std::move(json)};
    if (!res.is_success()) return; // keep error result unchanged
    using EmissionList = std::vector<pipeline::dynamic::Emission>;
    auto list = res.has_value() ? res.move_value<EmissionList>() : EmissionList{}; // Contacts returns EmissionList on success
    list.push_back(std::move(e));
    res = ComputationResult(std::move(list));
  }
};

} // namespace lahuta::analysis::contacts
