#ifndef LAHUTA_PIPELINE_DSL_TAGS_HPP
#define LAHUTA_PIPELINE_DSL_TAGS_HPP

namespace lahuta::pipeline::dsl {

struct thread_safe_t   { static constexpr bool value = true;  };
struct thread_unsafe_t { static constexpr bool value = false; };
inline constexpr thread_safe_t   thread_safe{};
inline constexpr thread_unsafe_t thread_unsafe{};

} // namespace lahuta::pipeline::dsl

#endif // LAHUTA_PIPELINE_DSL_TAGS_HPP
