/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr std::size_t align = alignof(std::string); alignas(align) char buf[sizeof(std::string)];
 *   auto* p = new (buf) std::string{"besian"}; p->append("sejdiu").append("@gmail.com");
 *   std::string r = *p; p->~basic_string(); return r;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_DSL_FUNCTION_TRAITS_HPP
#define LAHUTA_PIPELINE_DSL_FUNCTION_TRAITS_HPP

#include <functional>
#include <tuple>

namespace util {

/// fallback: T => traits<decltype(&T::operator())>
template <typename T, typename = void>
struct function_traits : function_traits<decltype(&T::operator())> {};

/// plain fn: R(Args...) => result_type=R, arity=|Args|, arg_t<N>=Args[N]
template <typename R, typename... Args>
struct function_traits<R(Args...)> {
  using result_type = R;

  static constexpr std::size_t arity = sizeof...(Args);

  template <std::size_t N>
  struct arg {
    static_assert(N < arity, "util::function_traits: parameter index out of range");
    using type = typename std::tuple_element<N, std::tuple<Args...>>::type;
  };

  template <std::size_t N>
  using arg_t = typename arg<N>::type;
};

// clang-format off
/// fn ptr/ref: R(*)(Args...)/R(&)(Args...) => <R(Args...)>
template <typename R, typename... Args> struct function_traits<R (*)(Args...)> : function_traits<R(Args...)> {};
template <typename R, typename... Args> struct function_traits<R (&)(Args...)> : function_traits<R(Args...)> {};

/// cv-qualified fn ptr: R(*const/volatile/const volatile)(Args...) => <R(Args...)>
template <typename R, typename... Args> struct function_traits<R (*const)(Args...)>          : function_traits<R(Args...)> {};
template <typename R, typename... Args> struct function_traits<R (*volatile)(Args...)>       : function_traits<R(Args...)> {};
template <typename R, typename... Args> struct function_traits<R (*const volatile)(Args...)> : function_traits<R(Args...)> {};
// clang-format on

/// std::function: std::function<R(Args...)> => <R(Args...)>
template <typename R, typename... Args>
struct function_traits<std::function<R(Args...)>> : function_traits<R(Args...)> {};

/// member-fn ptr: R(C::*)(Args...)[cv][ref][noexcept] => <R(Args...)>
#define DEF_TRAITS_CVREF_NOEX(fn_qual)                                                                       \
  template <typename C, typename R, typename... Args>                                                        \
  struct function_traits<R (C::*)(Args...) fn_qual> : function_traits<R(Args...)> {};

DEF_TRAITS_CVREF_NOEX()
DEF_TRAITS_CVREF_NOEX(const)
DEF_TRAITS_CVREF_NOEX(volatile)
DEF_TRAITS_CVREF_NOEX(const volatile)
DEF_TRAITS_CVREF_NOEX(&)
DEF_TRAITS_CVREF_NOEX(const &)
DEF_TRAITS_CVREF_NOEX(volatile &)
DEF_TRAITS_CVREF_NOEX(const volatile &)
DEF_TRAITS_CVREF_NOEX(&&)
DEF_TRAITS_CVREF_NOEX(const &&)
DEF_TRAITS_CVREF_NOEX(volatile &&)
DEF_TRAITS_CVREF_NOEX(const volatile &&)
DEF_TRAITS_CVREF_NOEX(noexcept)
DEF_TRAITS_CVREF_NOEX(const noexcept)
DEF_TRAITS_CVREF_NOEX(volatile noexcept)
DEF_TRAITS_CVREF_NOEX(const volatile noexcept)
DEF_TRAITS_CVREF_NOEX(& noexcept)
DEF_TRAITS_CVREF_NOEX(const & noexcept)
DEF_TRAITS_CVREF_NOEX(volatile & noexcept)
DEF_TRAITS_CVREF_NOEX(const volatile & noexcept)
DEF_TRAITS_CVREF_NOEX(&& noexcept)
DEF_TRAITS_CVREF_NOEX(const && noexcept)
DEF_TRAITS_CVREF_NOEX(volatile && noexcept)
DEF_TRAITS_CVREF_NOEX(const volatile && noexcept)

#undef DEF_TRAITS_CVREF_NOEX

/// alias: function_arg_t<F,N> = function_traits<F>::arg_t<N>
template <typename F, std::size_t N>
using function_arg_t = typename function_traits<F>::template arg_t<N>;

} // namespace util

#endif // LAHUTA_PIPELINE_DSL_FUNCTION_TRAITS_HPP
