/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string s;
 *   for (auto& p : parts) s += *std::addressof(p);
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_ENTITIES_SEARCH_CONFIG_HPP
#define LAHUTA_ENTITIES_SEARCH_CONFIG_HPP

#include <cstddef>

#include "entities/interaction_types.hpp"

// clang-format off
namespace lahuta::search {

enum class SearchAlgorithm { Brute, Grid };

//
// We store distance_max, but for robustness we overwrite this value using the contact-specific
// parameters. It's also not clear what, if any, performance benefit we get from better-guessing
// our allocations. Even more unclear is if this approach is worth it compared to just running
// the predicate-check for loop twice: once for count, and once for the actual checks.
//
struct SearchOptions {
  double distance_max        = 4.0f;
  float hit_reserve_factor   = 0.1f;
  float sel_reserve_factor_a = 0.2f;
  float sel_reserve_factor_b = sel_reserve_factor_a;
};

constexpr SearchOptions make_search_opts(Category c) {
  switch (c) {
    case Category::HydrogenBond:
      return { 4.1, 0.7f, 0.2f, 0.2f };
    case Category::WeakHydrogenBond:
      return { 4.1, 0.1f, 0.2f, 0.2f };
    case Category::Hydrophobic:
      return { 4.0, 0.5f, 0.3f, 0.3f };
    case Category::Halogen:
      return { 6.0, 0.0f, 0.0f, 0.3f };
    case Category::CationPi:
      return { 6.0, 0.1f, 0.5f, 1.0f };
    case Category::PiStacking:
      return { 6.0, 0.1f, 1.0f, 1.0f };
    case Category::MetalCoordination:
      return { 3.0, 0.0f, 0.0f, 0.3f };
    default:
      return {};
  }
}

} // namespace lahuta::search

#endif // LAHUTA_ENTITIES_SEARCH_CONFIG_HPP
