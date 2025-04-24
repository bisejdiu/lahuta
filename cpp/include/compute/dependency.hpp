#pragma once

#include "compute/label.hpp"
#include "compute/result.hpp"
#include <tuple>
#include <type_traits>
#include <vector>

// clang-format off
namespace lahuta::topology::compute {

struct NoDependencies {
  static std::vector<ComputationLabel> labels() { return {}; }
};

template <typename ComputationT, typename ResultT>
struct Dependency {
  using computation_type = ComputationT;
  using result_type = ResultT;

  static ComputationLabel label() { return ComputationT::label; }

  /// extract typed result
  static ResultT result(const ComputationResult& r) {
    return r.get_value<ResultT>();
  }
};

template <typename... Deps>
struct Dependencies {
  static std::vector<ComputationLabel> labels() {
    return {Deps::label()...};
  }
};

template <typename ParamT>
struct ParameterDependencies {
  using type = std::tuple<>;
};

template <typename ParamT>
using parameter_dependencies_t = typename ParameterDependencies<ParamT>::type;

/// Represents a single dependency relationship
template <typename ComputationT, typename ResultT>
struct DependencyInfo {
  using computation_type = ComputationT;
  using result_type = ResultT;

  static ComputationLabel label() { return ComputationT::label; }

  template <typename ParamT> static void update(ParamT &param, const ComputationResult &result) {
    if constexpr (std::is_same_v<void, ResultT>) {} 
    else { param.template update_from<ComputationT>(result.get_value<ResultT>()); }
  }
};

/// create a vector of ComputationLabels from a tuple of DependencyInfo
template <typename Tuple, std::size_t... Is>
std::vector<ComputationLabel> dependency_labels_impl(std::index_sequence<Is...>) {
  return {std::tuple_element_t<Is, Tuple>::label()...};
}

template <typename Tuple>
std::vector<ComputationLabel> dependency_labels() {
  return dependency_labels_impl<Tuple>(std::make_index_sequence<std::tuple_size_v<Tuple>>{});
}

/// Mixin class to add dependency handling to parameters
template <typename Derived>
class WithDependencies {
public:
  // Get dependency labels for this parameter type
  static std::vector<ComputationLabel> get_dependency_labels() {
    return dependency_labels<parameter_dependencies_t<Derived>>();
  }

  // Process a computation result
  static void process_dependency(Derived &param, const ComputationLabel &label, const ComputationResult &result) {
    process_impl<Derived>(
        param,
        label,
        result,
        std::make_index_sequence<std::tuple_size_v<parameter_dependencies_t<Derived>>>());
  }

private:
  // Match dependency by label and invoke the correct handler
  template <typename T, std::size_t... Is>
  static void process_impl(
      T &param, const ComputationLabel &label, const ComputationResult &result, std::index_sequence<Is...>) {
    // Fold expression to check each dependency
    (try_process_one<T, std::tuple_element_t<Is, parameter_dependencies_t<T>>>(param, label, result) || ...);
  }

  // Try to process a single dependency if labels match
  template <typename T, typename DepInfo>
  static bool try_process_one(T &param, const ComputationLabel &label, const ComputationResult &result) {
    if (label == DepInfo::label()) {
      DepInfo::template update<T>(param, result);
      return true;
    }
    return false;
  }
};

} // namespace lahuta::topology::compute
