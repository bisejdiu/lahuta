#ifndef LAHUTA_COMPUTE_DEPENDENCY_HPP
#define LAHUTA_COMPUTE_DEPENDENCY_HPP

#include <vector>

#include "compute/label.hpp"
#include "compute/result.hpp"

namespace lahuta::compute {

struct UnitComputation {
  static std::vector<ComputationLabel> labels() { return {}; }
};

template <typename ComputationT, typename ResultT>
struct Dependency {
  using computation_type = ComputationT;
  using result_type      = ResultT;

  static ComputationLabel label() { return ComputationT::label; }

  /// extract typed result
  static ResultT result(const ComputationResult &r) { //
    return r.get_value<ResultT>();
  }
};

template <typename... Deps>
struct Dependencies {
  static std::vector<ComputationLabel> labels() { //
    return {Deps::label()...};
  }
};

} // namespace lahuta::compute

#endif // LAHUTA_COMPUTE_DEPENDENCY_HPP
