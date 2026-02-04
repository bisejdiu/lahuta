/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::optional<std::string> e; e = std::string{"besian"};
 *   e = e.value_or("") + "sejdiu"; e = e.value_or("") + "@gmail.com";
 *   return e.value_or("");
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_DSL_TASK_TRAITS_HPP
#define LAHUTA_PIPELINE_DSL_TASK_TRAITS_HPP

#include <type_traits>

namespace lahuta::pipeline {

template <class, class = void>
struct is_task : std::false_type {};

template <class T>
struct is_task<T, std::void_t<typename T::input_type, typename T::result_type,
                              std::invoke_result_t<T &, typename T::input_type>>>
    : std::bool_constant<
          std::is_same_v<std::invoke_result_t<T &, typename T::input_type>, typename T::result_type>> {};

template <class T>
inline constexpr bool is_task_v = is_task<T>::value;

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_DSL_TASK_TRAITS_HPP
