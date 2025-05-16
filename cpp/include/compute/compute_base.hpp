#pragma once

#include "context.hpp"
#include "label.hpp"
#include "parameters.hpp"
#include "result.hpp"
#include <cstddef>
#include <memory>
#include <type_traits>
#include <vector>

namespace lahuta::topology::compute {

namespace detail {

/// trait checks if a computation type has dependencies defined
template <typename T, typename = void>
struct has_dependencies : std::false_type {};

template <typename T>
struct has_dependencies<T, std::void_t<typename T::dependencies>> : std::true_type {};

/// trait gets the dependencies of a computation type
template <typename T, typename = void>
struct dependencies_of {
  static std::vector<ComputationLabel> labels() { return {}; }
};

template <typename T>
struct dependencies_of<T, std::enable_if_t<has_dependencies<T>::value>> {
  static std::vector<ComputationLabel> labels() { return T::dependencies::labels(); }
};

} // namespace detail

template <typename DataT, Mut M>
class Computation {
public:
  virtual ~Computation() = default;

  /// Core execution method with base parameters
  virtual ComputationResult execute(DataContext<DataT, M>& context, const ParameterInterface& params) = 0;

  virtual std::unique_ptr<ParameterInterface>  get_parameters()   const = 0;
  virtual const ComputationLabel&              get_label()        const = 0;
  virtual std::vector<ComputationLabel>        get_dependencies() const = 0;

  // execution control
  bool is_enabled() const { return enabled_; }
  void set_enabled(bool enabled) { enabled_ = enabled; }

private:
  bool enabled_ = true;
};

} // namespace lahuta::topology::compute 
