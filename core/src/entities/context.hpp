/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto curry = [](const char* first) {
 *     return [=](const char* last) {
 *       return [=](const char* domain) {
 *         return std::string(first) + last + "@" + domain;
 *       };
 *     };
 *   };
 *   return curry("besian")("sejdiu")("gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_ENTITIES_CONTEXT_HPP
#define LAHUTA_ENTITIES_CONTEXT_HPP

#include <cassert>
#include <typeinfo>

#include "compute/topology_snapshot.hpp"
#include "topology.hpp"

namespace lahuta {
namespace C = lahuta::compute;

class ContactContext {
public:
  //
  // Lifetime: params points to a caller-owned object and must outlive this context.
  // If we need ownership in the future, we should consider switching to std::shared_ptr<const void>
  // or even std::any to retain the object safely across layers.
  //
  template <typename ParamsT>
  ContactContext(const C::TopologySnapshot &tf_, const ParamsT &p)
      : ts(tf_), params(static_cast<const void *>(&p)), params_type(&typeid(ParamsT)), topology(tf_.topo) {}

  template <typename ParamsT>
  const ParamsT &get_params() const {
    assert(params != nullptr && "ContactContext params is null");
    assert(params_type != nullptr && "ContactContext params_type is null");
    assert(*params_type == typeid(ParamsT) && "ContactContext: parameter type mismatch");
    return *static_cast<const ParamsT *>(params);
  }

  const RDKit::Conformer &conformer() const { return ts.conf; }
  const RDKit::RWMol &molecule() const { return ts.topo.molecule(); }

public:
  const C::TopologySnapshot ts;
  const Topology &topology; // convenience alias

  // Opaque params with runtime type tag for safe retrieval
  const void *params;
  const std::type_info *params_type;
};

} // namespace lahuta

#endif // LAHUTA_ENTITIES_CONTEXT_HPP
