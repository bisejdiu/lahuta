#ifndef LAHUTA_SEARCH_HPP
#define LAHUTA_SEARCH_HPP

#include <optional>
#include <type_traits>
#include <vector>

#include "entities.hpp"
/*#include "entity.hpp"*/
#include "nsgrid.hpp"

// clang-format off

namespace lahuta {

template <typename EV1, typename EV2>
using EnableIfEntityCollections =
    std::enable_if_t<
        std::is_base_of_v<EntityCollection<typename EV1::value_type>, EV1> &&
        std::is_base_of_v<EntityCollection<typename EV2::value_type>, EV2>>;

template <typename EV>
using EnableIfEntityCollection =
    std::enable_if_t<
        std::is_base_of_v<EntityCollection<typename EV::value_type>, EV>>;

//
// Using std::optional seems like a good idea, but it cannot hold references, and would force us to use
// a reference wrapper. Further, while it would simplify some of the code (constructors) it would
// obfuscate the intend of the code in distinguishing between self (enforce i < j) and non-self searches.
//
template <typename EV1, typename EV2, typename = EnableIfEntityCollections<EV1, EV2>>
class BruteForce {
public:
  BruteForce(const EV1 &ev_)                   : ev1(ev_), ev2(ev_), is_self_search(true) {}
  BruteForce(const EV1 &ev1_, const EV2 &ev2_) : ev1(ev1_), ev2(ev2_) {}

  NSResults search(double radius) { return search_entities(radius); }

private:
  NSResults search_entities(double radius) {

    NSResults results;
    double dist_sq = radius * radius;

    for (size_t idx_a = 0; idx_a < ev1.size(); ++idx_a) {
      for (size_t idx_b = 0; idx_b < ev2.size(); ++idx_b) {

        if (is_self_search && !is_valid_pair(idx_b, idx_a)) {
          continue;
        }

        auto dist_result = get_dist_within_cutoff(ev1[idx_a], ev2[idx_b], dist_sq);
        if (dist_result.has_value()) {
          results.add(idx_a, idx_b, dist_result.value());
        }
      }
    }

    return results;
  }

  template <typename EntityA, typename EntityB>
  std::optional<double> get_dist_within_cutoff(const EntityA &a, const EntityB &b, double dist_sq) {

    auto pos_as = get_pos(a);
    auto pos_bs = get_pos(b);

    for (const auto pos_a : pos_as) {
      for (const auto pos_b : pos_bs) {
        double current_dist_sq = ((*pos_a) - (*pos_b)).lengthSq();
        if (current_dist_sq < dist_sq) {
          return current_dist_sq;
        }
      }
    }
    return std::nullopt;
  }

  bool is_valid_pair(int i, int j) {
    if (i >= j) return false;
    return true;
  };

  inline static auto get_pos(const Entity& entity) -> std::vector<const RDGeom::Point3D*> {
    return { &entity.get_center() };
  }

private:
  const EV1 &ev1;
  const EV2 &ev2;
  bool is_self_search = false;
};

template <typename EV1, typename EV2, typename = EnableIfEntityCollections<EV1, EV2>>
struct SearchStrategy {
  static constexpr bool prefer_grid(const EV1 &ev1, const EV2 &ev2, double) {
    return ev1.get_data().size() * ev2.get_data().size() > 500;
  }

  static constexpr bool size_is_valid(const EV1 &ev1, const EV2 &ev2) {
    return ev1.get_data().size() > 0 && ev2.get_data().size() > 0;
  }
};

template <typename EV1, typename EV2, typename = EnableIfEntityCollection<EV1>>
class SearchImpl {
public:
  NSResults brute_force(const EV1 &ev1, double radius) const {
    BruteForce<EV1, EV1> brute_force(ev1);
    return brute_force.search(radius);
  }

  NSResults brute_force(const EV1 &ev1, const EV2 &ev2, double radius) const {
    BruteForce<EV1, EV2> brute_force(ev1, ev2);
    return brute_force.search(radius);
  }

  NSResults _grid_impl(const EV1 &ev1, const EV2 &ev2, double radius) const {
    auto grid = FastNS(ev2.positions());
    auto ok = grid.build(radius);
    return ok ? grid.search(ev1.positions()) : brute_force(ev1, ev2, radius);
  }

  NSResults _grid_impl(const EV1 &ev1, double radius) const {
    auto grid = FastNS(ev1.positions());
    auto ok = grid.build(radius);
    return ok ? grid.self_search() : brute_force(ev1, radius);
  }

  NSResults grid(const EV1 &ev1, const std::optional<EV2> &ev2, double radius) const {
    if (radius <= 0.0) throw std::invalid_argument("Radius must be greater than 0.0");

    if (ev2.has_value())
      return _grid_impl(ev1, ev2.value(), radius);

    return _grid_impl(ev1, radius);
  }
};

class EntityNeighborSearch {
public:
template <typename EV1, typename EV2, typename = EnableIfEntityCollections<EV1, EV2>>
  static NSResults search(const EV1 &ev1, const EV2 &ev2, double radius) {
    if (!SearchStrategy<EV1, EV2>::size_is_valid(ev1, ev2)) return NSResults();
    SearchImpl<EV1, EV2> impl;
    return SearchStrategy<EV1, EV2>::prefer_grid(ev1, ev2, radius)
               ? impl.grid(ev1, ev2, radius)
               : impl.brute_force(ev1, ev2, radius);
  }

  template <typename EV, typename = EnableIfEntityCollection<EV>>
  static NSResults search(const EV &ev, double radius) {
    if (!SearchStrategy<EV, EV>::size_is_valid(ev, ev)) return NSResults();
    SearchImpl<EV, EV> impl;
    return SearchStrategy<EV, EV>::prefer_grid(ev, ev, radius)
               ? impl.grid(ev, std::nullopt, radius)
               : impl.brute_force(ev, radius);
  }
};

} // namespace lahuta

#endif // LAHUTA_SEARCH_HPP
