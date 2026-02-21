/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr std::size_t align = alignof(std::string); alignas(align) char buf[sizeof(std::string)];
 *   auto* p = new (buf) std::string{"besian"}; p->append("sejdiu").append("@gmail.com");
 *   std::string r = *p; p->~basic_string(); return r;
 * }();
 *
 */

#ifndef LAHUTA_ENTITIES_SEARCH_HPP
#define LAHUTA_ENTITIES_SEARCH_HPP

#include <cstdint>
#include <vector>

#include <rdkit/Geometry/point.h>

#include "entities/search/config.hpp"
#include "entities/search/hit_buffer.hpp"
#include "entities/search/provider.hpp"
#include "entities/search/search_impl.hpp"
#include "logging/logging.hpp"
#include "utils/span.hpp"

// clang-format off
namespace lahuta { struct ContactContext;}
namespace lahuta::search {

static constexpr std::size_t GRID_THRESHOLD = 5000;

struct NoTester {
  template <typename... Args>
  constexpr bool operator()(Args &&...) const noexcept { return true; }
};

inline constexpr SearchAlgorithm chose_algo_heuristic(std::size_t size_a, std::size_t size_b) noexcept {
  const bool use_grid = size_a * size_b > GRID_THRESHOLD; // trying to be smart here just makes it easy to shoot yourself in the foot
  return use_grid ? SearchAlgorithm::Grid : SearchAlgorithm::Brute;
}

struct NeighborResult {
  std::vector<Hit> hits;
  std::vector<uint32_t> idxs_a;
  std::vector<uint32_t> idxs_b;
  bool is_self_search;
};

// Materialize a tested set of entity indices
template <typename RecordT, typename PredT>
std::vector<uint32_t> materialize_selection(const std::vector<RecordT> &records, const PredT &pred, const float reserve_factor) {

  std::vector<uint32_t> result;
  result.reserve(records.size() * reserve_factor);

  for (uint32_t i = 0; i < records.size(); ++i) {
    if (pred(records[i])) result.push_back(i);
  }

  return result;
}

template<
    bool SelfSearch,
    typename CoordT,
    typename RecordT1, typename PredT1,
    typename RecordT2, typename PredT2,
    typename Tester
>
NeighborResult neighbour_search(
  const CoordT &coord_provider,
  const std::vector<RecordT1> &recs_a, const PredT1 &pred_a,
  const std::vector<RecordT2> &recs_b, const PredT2 &pred_b,
  const search::SearchOptions opts, Tester &&tester, const ContactContext& ctx) {

  static_assert(is_coord_provider_v<CoordT>, "CoordProvider must expose get_a()/get_b() to Point3D&");

  NeighborResult result;
  result.is_self_search = SelfSearch;

  result.idxs_a = materialize_selection(recs_a, pred_a, opts.sel_reserve_factor_a);
  if constexpr (SelfSearch) {
    result.idxs_b = result.idxs_a;
  } else {
    result.idxs_b = materialize_selection(recs_b, pred_b, opts.sel_reserve_factor_b);
  }

  Logger::get_logger()->debug("NeighborSearch: {} search, {} entities A, {} entities B, max_distance={:.2f}",
                              SelfSearch ? "Self" : "Cross",
                              result.idxs_a.size(), result.idxs_b.size(), opts.distance_max);

  if (result.idxs_a.empty() || result.idxs_b.empty()) return result;

  span<uint32_t> span_a(result.idxs_a);
  span<uint32_t> span_b(result.idxs_b);

  HitBuffer buffer;
  buffer.reserve(std::min(result.idxs_a.size(), result.idxs_b.size()) * opts.hit_reserve_factor);

  [&]() {
    float radius_sq = opts.distance_max * opts.distance_max;
    SearchAlgorithm algo = chose_algo_heuristic(span_a.size(), span_b.size());

    Logger::get_logger()->debug("NeighborSearch: Using {} search algorithm for {}x{} entities",
                                algo == SearchAlgorithm::Grid ? "Grid" : "Brute Force",
                                span_a.size(), span_b.size());

    if (algo == SearchAlgorithm::Grid) {

      auto &grid_span = SelfSearch ? span_a : span_b;
      GridSearch<SelfSearch, CoordT> grid{coord_provider, grid_span, opts.distance_max};

      if (grid(span_a, span_b, buffer, std::forward<Tester>(tester), ctx)) return; // success
      Logger::get_logger()->warn("Grid search failed, falling back to brute force");
    }
    BruteForceSearch<SelfSearch, CoordT> brute{coord_provider, radius_sq};
    brute(span_a, span_b, buffer, std::forward<Tester>(tester), ctx);
  }();

  result.hits = buffer.release();

  Logger::get_logger()->debug("NeighborSearch: Found {} potential contacts", result.hits.size());

  return result;
}

} // namespace lahuta::search

#endif // LAHUTA_ENTITIES_SEARCH_HPP
