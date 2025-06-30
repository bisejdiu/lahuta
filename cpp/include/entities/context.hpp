#ifndef LAHUTA_ENTITIES_CONTEXT_HPP
#define LAHUTA_ENTITIES_CONTEXT_HPP

#include "topology.hpp"

namespace lahuta {

struct ContactContext {
  const Topology& topology;
  const void* params; // not happy with this type erasure.

  template<typename ParamsT>
  ContactContext(const Topology& topo, const ParamsT& p) : topology(topo), params(&p) {}

  template<typename ParamsT>
  const ParamsT& get_params() const {
    return *static_cast<const ParamsT*>(params);
  }

  const RDKit::RWMol& molecule() const { return topology.molecule(); }
};

} // namespace lahuta

#endif // LAHUTA_ENTITIES_CONTEXT_HPP
