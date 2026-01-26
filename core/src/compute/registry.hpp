#ifndef LAHUTA_COMPUTE_REGISTRY_HPP
#define LAHUTA_COMPUTE_REGISTRY_HPP

#include <algorithm>
#include <stdexcept>
#include <utility>

#include "compute/compute_base.hpp"
#include "node.hpp"

namespace lahuta::compute {

template <typename D, Mut M = Mut::ReadWrite>
class Registry {
public:
  /// add a computation to the registry, returning its index
  int add(std::unique_ptr<Computation<D, M>> computation) {
    if (size_ >= MAX_N_COMPUTATIONS) throw std::runtime_error("MAX_N_COMPUTATIONS reached");
    int node_idx = find(computation->get_label());
    if (node_idx >= 0) throw std::runtime_error("duplicate computation");

    node_idx   = size_++;
    auto &node = nodes[node_idx];
    node.tag   = computation->get_label();
    node.deps  = 0; // filled later
    node.impl  = std::move(computation);
    node.proto = node.impl->get_parameters();
    return node_idx;
  }

  /// seal the registry, i.e. fill in the (forward & reverse) dependencies
  void seal() {
    for (int node_idx = 0; node_idx < size_; ++node_idx)
      for (auto dep_label : nodes[node_idx].impl->get_dependencies()) {
        int dep_idx = find(dep_label);
        if (dep_idx < 0) throw std::runtime_error("missing dependency");

        nodes[node_idx].deps |= (Mask{1} << dep_idx);
        nodes[dep_idx].rdeps |= (Mask{1} << node_idx);
      }
  }

  /// find a computation by its label, returning its index or -1 if not found
  int find(ComputationLabel label) const {
    for (int node_idx = 0; node_idx < size_; ++node_idx)
      if (nodes[node_idx].tag == label) return node_idx;
    return -1;
  }

  ComputeNode<D, M> &operator[](int i) { return nodes[i]; }
  const ComputeNode<D, M> &operator[](int i) const { return nodes[i]; }
  int size() const { return size_; }

private:
  std::array<ComputeNode<D, M>, MAX_N_COMPUTATIONS> nodes{};
  u8 size_ = 0;
};

} // namespace lahuta::compute

#endif // LAHUTA_COMPUTE_REGISTRY_HPP
