#ifndef LAHUTA_PIPELINE_TASK_TRAITS_HPP
#define LAHUTA_PIPELINE_TASK_TRAITS_HPP

#include <type_traits>
#include <utility>

// clang-format off
namespace lahuta::pipeline::dsl {

template <typename, typename = void>
struct is_task : std::false_type {};

template <typename T>
struct is_task<
  T, std::void_t<
       typename T::input_type, typename T::result_type,
       decltype(std::declval<T>()(std::declval<typename T::input_type>()))>>
  : std::bool_constant<std::is_same_v<
       decltype(std::declval<T>()(std::declval<typename T::input_type>())), typename T::result_type>> {};

template <typename T>
inline constexpr bool is_task_v = is_task<T>::value;

} // namespace lahuta::pipeline::dsl

#endif // LAHUTA_PIPELINE_TASK_TRAITS_HPP
