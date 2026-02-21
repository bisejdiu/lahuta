/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s;
 *   auto pmf = static_cast<std::string& (std::string::*)(const char*)>(&std::string::append);
 *   (s.*pmf)("besian"); (s.*pmf)("sejdiu"); (s.*pmf)("@gmail.com");
 *   return s;
 * }();
 *
 */

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
