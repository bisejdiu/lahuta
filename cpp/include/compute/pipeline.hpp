#pragma once

#include "node.hpp"
#include "registry.hpp"

namespace lahuta::topology::compute {

/// Represents a sequence of computations that respects dependency ordering.
/// Validates the execution order at construction time to ensure all dependencies
/// are satisfied before their dependent computations.
template <typename D, Mut M>
class Pipeline {
  ExecOrder plan_; // immutable
public:
  explicit Pipeline(const Registry<D, M> &registry, std::initializer_list<ComputationLabel> ordered) {
    // validate dependenxies in O(N) time
    Mask satisfied_deps = 0;
    for (auto label : ordered) {
      int node_idx = registry.find(label);
      if (node_idx < 0) throw std::runtime_error("unknown label in pipeline");

      // every dependency of this node must already be in 'satisfied_deps'
      Mask node_deps = registry[node_idx].deps;
      for (int dep_idx = 0; dep_idx < registry.size(); ++dep_idx) {
        if (node_deps & (1 << dep_idx) && !(satisfied_deps & (1 << dep_idx))) {
          throw std::runtime_error("pipeline order violates dependency for " + std::string(registry[node_idx].tag.to_string_view()));
        }
      }

      plan_.node_indices[plan_.size++] = static_cast<u8>(node_idx);
      satisfied_deps |= (1 << node_idx);
    }
  }
  const ExecOrder &plan() const noexcept { return plan_; }
};

} // namespace lahuta::topology::compute
