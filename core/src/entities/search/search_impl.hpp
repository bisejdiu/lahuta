#ifndef LAHUTA_ENTITIES_SEARCH_IMPL_HPP
#define LAHUTA_ENTITIES_SEARCH_IMPL_HPP

#include <vector>

#include "distances/neighbors.hpp"
#include "entities/context.hpp"
#include "entities/search/hit_buffer.hpp"
#include "entities/search/provider.hpp"
#include "spatial/nsresults.hpp"

// clang-format off
namespace lahuta::search {

template <bool SelfSearch, typename CoordProvider, typename Buffer = search::HitBuffer>
struct BruteForceSearch {
  static_assert(is_coord_provider_v<CoordProvider>, "CoordProvider must expose get_a()/get_b() to Point3D&");
  static_assert(is_hit_buffer_v<Buffer>, "Buffer must support push(i, j, d2, tester)");

  const CoordProvider &coord;
  float radius_sq;

  template <typename SpanA, typename SpanB, typename Tester>
  bool operator()(const SpanA &idxs_a, const SpanB &idxs_b, Buffer &hits, Tester &&f, const ContactContext& ctx) const {
    for (std::size_t ix_a = 0; ix_a < idxs_a.size(); ++ix_a) {
      const auto &a_pos = coord.get_a(idxs_a[ix_a]);

      std::size_t start_b = SelfSearch ? ix_a + 1 : 0;
      for (std::size_t ix_b = start_b; ix_b < idxs_b.size(); ++ix_b) {
        const auto &b_pos = SelfSearch ? coord.get_a(idxs_b[ix_b]) : coord.get_b(idxs_b[ix_b]);

        float dx = a_pos.x - b_pos.x, dy = a_pos.y - b_pos.y, dz = a_pos.z - b_pos.z;
        float d2 = dx * dx + dy * dy + dz * dz;
        if (d2 <= radius_sq) {
          hits.push(idxs_a[ix_a], SelfSearch ? idxs_a[ix_b] : idxs_b[ix_b], d2, std::forward<Tester>(f), ctx);
        }
      }
    }

    return true;
  }
};

template <bool SelfSearch, typename CoordProvider, typename Buffer = HitBuffer>
struct GridSearch {
  static_assert(is_coord_provider_v<CoordProvider>, "CoordProvider must expose get_a()/get_b() to Point3D&");
  static_assert(is_hit_buffer_v<Buffer>);

  const CoordProvider &coord;
  double radius;
  std::vector<RDGeom::Point3D> build_pts;

  template <typename SpanBuild>
  GridSearch(const CoordProvider &cp, const SpanBuild &build_idxs, double r) : coord(cp), radius(r) {
    build_pts.reserve(build_idxs.size());
    for (auto idx : build_idxs) {
      build_pts.push_back(SelfSearch ? coord.get_a(idx) : coord.get_b(idx));
    }
  }

  template <typename SpanA, typename SpanB, typename Tester>
  bool operator()(const SpanA &idxs_a, const SpanB &idxs_b, Buffer &hits, Tester &&f, const ContactContext& ctx) const {
    if constexpr (SelfSearch)
      return self_search(idxs_a, hits, std::forward<Tester>(f), ctx);

    return cross_search(idxs_a, idxs_b, hits, std::forward<Tester>(f), ctx);
  }

private:
  template <typename Span, typename Tester>
  bool self_search(const Span &idxs, Buffer &hits, Tester &&f, const ContactContext& ctx) const {

    dist::NeighborSearchOptions opts;
    opts.cutoff = radius;
    NSResults results = dist::neighbors_within_radius_self(build_pts, opts);

    for (size_t res = 0; res < results.size(); ++res) {
      auto [i, j] = results.get_pairs()[res];
      float d2 = results.get_distances()[res];
      if (i < idxs.size() && j < idxs.size()) {
        hits.push(idxs[i], idxs[j], d2, std::forward<Tester>(f), ctx);
      }
    }

    return true;
  }

  template <typename SpanA, typename SpanB, typename Tester> // A is the query
  bool cross_search(const SpanA &idxs_a, const SpanB &idxs_b, Buffer &hits, Tester &&f, const ContactContext& ctx) const {

    std::vector<RDGeom::Point3D> q_points;
    q_points.reserve(idxs_a.size());
    for (auto idx : idxs_a) {
      q_points.push_back(coord.get_a(idx));
    }

    dist::NeighborSearchOptions opts;
    opts.cutoff = radius;
    NSResults results = dist::neighbors_within_radius_cross_fastns(q_points, build_pts, opts);

    for (size_t res = 0; res < results.size(); ++res) {
      auto [qi, bi] = results.get_pairs()[res]; // query, build
      float d2 = results.get_distances()[res];
      if (qi < idxs_a.size() && bi < idxs_b.size()) {
        hits.push(idxs_a[qi], idxs_b[bi], d2, std::forward<Tester>(f), ctx);
      }
    }

    return true;
  }
};

} // namespace lahuta::search

#endif // LAHUTA_ENTITIES_SEARCH_IMPL_HPP
