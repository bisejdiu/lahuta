#ifndef LAHUTA_SEARCH_HPP
#define LAHUTA_SEARCH_HPP

#include "Geometry/point.h"
#include "GraphMol/Conformer.h"
#include "contacts/contacts.hpp"
#include "contacts/features.hpp"
#include "nsgrid.hpp"
#include "rings.hpp"
#include <vector>

namespace lahuta {

/*inline std::vector<const RDKit::Atom *> get_atoms(const AtomData &a) {*/
/*  return {a.atom};*/
/*}*/
/**/
/*inline std::vector<const RDKit::Atom *> get_atoms(const Feature &f) {*/
/*  return f.members;*/
/*}*/
/**/
/*inline std::vector<const RDKit::Atom *> get_atoms(const RingData &f) {*/
/*  return f.atoms;*/
/*}*/

inline std::vector<const RDGeom::Point3D *> get_pos(const RDKit::Conformer &conf, const AtomData &a) {
  return {a.pos};
}

inline std::vector<const RDGeom::Point3D *> get_pos(const RDKit::Conformer &conf, const Feature &f) {
  auto pos = new RDGeom::Point3D(f.center);
  return {pos};
}

inline std::vector<RDGeom::Point3D *> get_pos(const RDKit::Conformer &conf, const RingData &f) {
  auto pos = new RDGeom::Point3D(f.center);
  return {pos};
}

template <typename EV1, typename EV2> // entity vector
class BruteForce {
public:
  BruteForce(const RDKit::Conformer &conf_, const EV1 &ev1_, const EV2 &ev2_)
      : conf(conf_), ev1(ev1_), ev2(ev2_) {
  }
  BruteForce(const RDKit::Conformer &conf_, const EV1 &ev_)
      : conf(conf_), ev1(ev_), ev2(ev_), is_self_search(true) {}

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

    auto pos_as = get_pos(conf, a);
    auto pos_bs = get_pos(conf, b);

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

private:
  const RDKit::Conformer &conf;
  const EV1 &ev1;
  const EV2 &ev2;
  bool is_self_search = false;
};

template <typename EV1, typename EV2> //
struct SearchStrategy {
  static constexpr bool prefer_grid(const EV1 &ev1, const EV2 &ev2, double) {
    return ev1.get_data().size() * ev2.get_data().size() > 1000;
  }
};

template <typename EV1, typename EV2> //
class SearchImpl {
  const RDKit::Conformer &conf;

public:
  explicit SearchImpl(const RDKit::Conformer &conf_) : conf(conf_) {}

  NSResults brute_force(const EV1 &ev1, const EV2 &ev2, double radius) const {
    BruteForce<EV1, EV2> brute_force(conf, ev1, ev2);
    return brute_force.search(radius);
  }

  NSResults grid(const EV1 &ev1, const EV2 &ev2, double radius) const {
    auto grid = FastNS::create(ev2.positions(), radius);
    return grid.is_valid() ? grid.search(ev1.positions()) : brute_force(ev1, ev2, radius);
  }

  NSResults brute_force(const EV1 &ev1, double radius) const {
    BruteForce<EV1, EV1> brute_force(conf, ev1);
    return brute_force.search(radius);
  }

  NSResults grid(const EV1 &ev1, double radius) const {
    auto grid = FastNS::create(ev1.positions(), radius);
    return grid.is_valid() ? grid.self_search() : brute_force(ev1, radius);
  }
};

class EntityNeighborSearch {
public:
  explicit EntityNeighborSearch(const RDKit::Conformer &conf_) : conf(conf_) {}

  template <typename EV1, typename EV2> //
  NSResults search(const EV1 &ev1, const EV2 &ev2, double radius) {
    SearchImpl<EV1, EV2> impl(conf);
    return SearchStrategy<EV1, EV2>::prefer_grid(ev1, ev2, radius) //
               ? impl.grid(ev1, ev2, radius)
               : impl.brute_force(ev1, ev2, radius);
  }

  template <typename EV> NSResults search(const EV &ev, double radius) {
    SearchImpl<EV, EV> impl(conf);
    return SearchStrategy<EV, EV>::prefer_grid(ev, ev, radius) //
               ? impl.grid(ev, radius)
               : impl.brute_force(ev, radius);
  }

private:
  const RDKit::Conformer &conf;
};

} // namespace lahuta

#endif // LAHUTA_SEARCH_HPP
