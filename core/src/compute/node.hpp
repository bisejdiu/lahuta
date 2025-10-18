#ifndef LAHUTA_COMPUTE_NODE_HPP
#define LAHUTA_COMPUTE_NODE_HPP

#include "compute/_defs.hpp"
#include "compute/compute_base.hpp"
#include "compute/label.hpp"
#include "parameters.hpp"
#include "result.hpp"
#include <array>
#include <cassert>
#include <cstdint>
#include <memory>

// clang-format off
namespace lahuta::topology::compute {

constexpr std::size_t MAX_N_COMPUTATIONS = 64;
using Mask = std::uint64_t;

struct ExecOrder {
  std::array<u8, MAX_N_COMPUTATIONS> node_indices{};
  u8 size = 0;
};

template <typename DataT, Mut M>
struct ComputeNode {
  ComputationLabel tag{""};
  std::unique_ptr<Computation<DataT, M>> impl; // sole owner
  std::unique_ptr<ParameterInterface> proto;   // parameter prototype
  bool forced_disabled = false;      // if true, the computation is hard disabled
  Mask deps     = 0;                 // bit graph
  Mask rdeps    = 0;                 // computed once
  bool enabled  = true;              // decides if a computation participates
  // these two cache the last run
  bool done     = false;             // if the computation was run
  ComputationResult res;             // the result of the last run
  bool postprocessed = false;        // if on_complete was invoked on the result
};

} // namespace lahuta::topology::compute

#endif // LAHUTA_COMPUTE_NODE_HPP
