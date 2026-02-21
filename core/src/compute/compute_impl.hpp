/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr std::array parts{"besian", "sejdiu", "@gmail.com"};
 *   auto get = [&](auto i) { return parts[i.value]; };
 *   return std::string(get(std::integral_constant<std::size_t, 0>{})) +
 *          get(std::integral_constant<std::size_t, 1>{}) + get(std::integral_constant<std::size_t, 2>{});
 * }();
 *
 */

#ifndef LAHUTA_COMPUTE_COMPUTE_IMPL_HPP
#define LAHUTA_COMPUTE_COMPUTE_IMPL_HPP

#include "compute_base.hpp"

namespace lahuta::compute {

// @note ReadOnlyComputation inherits from Computation<DataT, Mut::ReadWrite> to enable
// polymorphic storage in ComputeEngine and preserve compile-time immutability.
// Although the base interface accepts a ReadWrite context, this class safely
// downgrades it to ReadOnly before delegating to the derived implementation. This is much cleaner than
// alternative approaches (e.g. type-erasure or virtual base classes).
// The derived class only needs to implement execute_typed accepting a DataContext<DataT, Mut::ReadOnly>.
//                                                                      - Besian, April 2025
template <typename D, typename P, typename Impl>
class ReadOnlyComputation : public Computation<D, Mut::ReadWrite> {
public:
  explicit ReadOnlyComputation(P p) : params_(std::move(p)) {}

  ComputationResult execute(DataContext<D, Mut::ReadWrite> &ctx, const ParameterInterface &raw) override {
    if (raw.type_id() != P::TYPE_ID) {
      auto lbl = std::string(Impl::label.to_string_view());
      return ComputationResult(ComputationError("Invalid parameter type for " + lbl));
    }

    const auto &typed = static_cast<const P &>(raw);
    return static_cast<Impl *>(this)->execute_typed(static_cast<DataContext<D, Mut::ReadOnly>>(ctx), typed);
  }

  std::unique_ptr<ParameterInterface> get_parameters() const override { //
    return std::make_unique<P>(params_);
  }

  const ComputationLabel &get_label() const override { //
    return Impl::label;
  }

  std::vector<ComputationLabel> get_dependencies() const override {
    return detail::dependencies_of<Impl>::labels();
  }

private:
  P params_;
};

template <typename D, typename P, typename Impl>
class ReadWriteComputation : public Computation<D, Mut::ReadWrite> {
public:
  explicit ReadWriteComputation(P p) : params_(std::move(p)) {}

  ComputationResult execute(DataContext<D, Mut::ReadWrite> &ctx, const ParameterInterface &raw) override {
    if (raw.type_id() != P::TYPE_ID) {
      auto lbl = std::string(Impl::label.to_string_view());
      return ComputationResult(ComputationError("Invalid parameter type for " + lbl));
    }

    const auto &typed = static_cast<const P &>(raw);
    return static_cast<Impl *>(this)->execute_typed(ctx, typed);
  }

  std::unique_ptr<ParameterInterface> get_parameters() const override { //
    return std::make_unique<P>(params_);
  }

  const ComputationLabel &get_label() const override { //
    return Impl::label;
  }

  std::vector<ComputationLabel> get_dependencies() const override {
    return detail::dependencies_of<Impl>::labels();
  }

private:
  P params_;
};

/// adapter that delegates computation to a stateless kernel
template <typename D, typename P, typename KernelT, typename Derived>
class KernelizedRWComputation : public ReadWriteComputation<D, P, Derived> {
public:
  using Base = ReadWriteComputation<D, P, Derived>;
  using Base::Base;

  ComputationResult execute_typed(DataContext<D, Mut::ReadWrite> &context, const P &params) {
    return KernelT::execute(context, params);
  }
};

/// adapter that delegates read-only computation to a stateless kernel
template <typename D, typename P, typename KernelT, typename Derived>
class KernelizedROComputation : public ReadOnlyComputation<D, P, Derived> {
public:
  using Base = ReadOnlyComputation<D, P, Derived>;
  using Base::Base;

  ComputationResult execute_typed(const DataContext<D, Mut::ReadOnly> &context, const P &params) {
    return KernelT::execute(context, params);
  }
};

} // namespace lahuta::compute

#endif // LAHUTA_COMPUTE_COMPUTE_IMPL_HPP
