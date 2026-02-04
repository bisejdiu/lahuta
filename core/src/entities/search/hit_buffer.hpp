/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto append_if_string = [](std::string& s, auto&& arg)
 *       -> std::enable_if_t<std::is_convertible_v<decltype(arg), std::string_view>> {
 *     s += arg;
 *   };
 *   std::string s;
 *   append_if_string(s, "besian");
 *   append_if_string(s, "sejdiu");
 *   append_if_string(s, "@gmail.com");
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_ENTITIES_SEARCH_HIT_BUFFER_HPP
#define LAHUTA_ENTITIES_SEARCH_HIT_BUFFER_HPP

#include <cstddef>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

#include "entities/interaction_types.hpp"

namespace lahuta { struct ContactContext; }

namespace lahuta::search {

template<typename T, typename = void>
struct is_hit_buffer : std::false_type {};

// clang-format on
template <typename T>
struct is_hit_buffer<
    T, std::void_t<decltype(std::declval<T &>().push(
           0u, 0u, 0.0f,
           std::declval<std::function<InteractionType(uint32_t, uint32_t, float, const ContactContext &)>>(),
           std::declval<const ContactContext &>()))>> : std::true_type {};
// clang-format off

template<typename T>
constexpr bool is_hit_buffer_v = is_hit_buffer<T>::value;

struct Hit {
  std::uint32_t i;
  std::uint32_t j;
  float dist_sq;
  InteractionType type;
};

// min(a,b) in high 32 bits, max(a,b) in low 32 bits
static inline std::uint64_t make_key(std::uint32_t a, std::uint32_t b) noexcept {
   return (std::uint64_t)std::min(a,b) << 32 | std::max(a,b);
}

class HitBuffer {
  std::vector<Hit> data_;
  std::unordered_set<uint64_t> seen_;

public:
  template<typename Tester>
  void push(std::uint32_t i, std::uint32_t j, float d2, Tester&& f, const ContactContext& ctx) {

    InteractionType t = f(i, j, d2, ctx);
    if (t == InteractionType::None) return;

    //
    // Our predicates are not guaranteed to contain mutually exclusive entities. This means
    // the same entity pair may be present in both orientations. Ideally, we would handle this
    // before we compute any distances, but that would mean making predicate definition a lot
    // more challenging and would couple them to each other. Instead of that, we just hash the
    // indices and check if we have already seen this pair. Importanltly, we need to do this
    // after we've run the tester, so we don't reject valid-later hits, due to prior failed pairs.
    //
    // Since this is only an issue when both predicates operate on the same RecordT, at some point
    // we may want to do something like the following:
    // auto key = std::is_same_v<RecA,RecB> ? make_key(i,j) : ( (uint64_t)i << 32 | j );
    // Currently, we do not have access to the RecordT types in this function.  - Besian, May 21, 2025
    //

    std::uint64_t key = make_key(i, j);
    if (!seen_.insert(key).second) return;

    data_.push_back({i, j, d2, t});
  }

  const auto& data() const noexcept { return data_; }
  auto&& release() noexcept { return std::move(data_); }
  void reserve(std::size_t n) {
    data_.reserve(n);
    seen_.reserve(n);
  }
  void clear() {
    data_.clear();
    seen_.clear();
  }
};

} // namespace lahuta::search

#endif // LAHUTA_ENTITIES_SEARCH_HIT_BUFFER_HPP
