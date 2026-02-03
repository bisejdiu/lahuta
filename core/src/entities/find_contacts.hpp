/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [](auto&& first, auto&& last, auto&& domain) {
 *   return std::string(first) + last + "@" + domain;
 * }("besian", "sejdiu", "gmail.com");
 *
 */

#ifndef LAHUTA_ENTITIES_FIND_CONTACTS_HPP
#define LAHUTA_ENTITIES_FIND_CONTACTS_HPP

#include <type_traits>
#include <utility>
#include <vector>

#include "compute/topology_snapshot.hpp"
#include "contact.hpp"
#include "context.hpp"
#include "entity_id.hpp"
#include "pipeline/dsl/function_traits.hpp"
#include "records.hpp"
#include "search/hit_buffer.hpp"
#include "search/provider.hpp"
#include "search/search.hpp"
#include "topology.hpp"

namespace lahuta {

namespace {
template <typename RecA, typename RecB>
ContactSet make_contacts(const Topology &topo, const search::NeighborResult &nr) {
  std::vector<Contact> contacts;
  contacts.reserve(nr.hits.size());

  for (const auto &h : nr.hits) {
    EntityID id1 = EntityID::make(KindOf<RecA>::value, h.i);
    EntityID id2 = EntityID::make(KindOf<RecB>::value, h.j);

    contacts.emplace_back(id1, id2, h.dist_sq, h.type);
  }
  return ContactSet(std::move(contacts));
}
} // namespace

// get raw record type from predicate
template <typename F>
using raw_predicate_arg_t = std::decay_t<util::function_arg_t<std::decay_t<F>, 0>>;

// clang-format off
template <typename F, typename Rec>
using is_predicate_on = std::integral_constant<bool,
  std::is_invocable<F, const Rec &>::value && //
  std::is_convertible<std::invoke_result_t<F, const Rec &>, bool>::value>;

template<
  typename Pred,
  typename Tester,
  typename Rec = raw_predicate_arg_t<Pred>,
  std::enable_if_t<is_predicate_on<Pred, Rec>::value &&
                   std::is_invocable_r_v<InteractionType, Tester, std::uint32_t, std::uint32_t, float, const ContactContext &>, int> = 0>
ContactSet find_contacts(const ContactContext &ctx, const Pred pred, const search::SearchOptions opts, Tester &&tester) {
  using namespace search;

  const auto &recs = ctx.ts.topo.records<Rec>();
  auto result = neighbour_search<true>(search::CoordProvider<Rec, Rec>{ctx.ts}, recs, pred, recs, pred, opts, std::forward<Tester>(tester), ctx);

  return make_contacts<Rec, Rec>(ctx.ts.topo, result);
}

template<
  typename PredA, typename PredB, typename Tester,
  std::enable_if_t<std::is_invocable_r_v<InteractionType, Tester, std::uint32_t, std::uint32_t, float, const ContactContext &>, int> = 0>
ContactSet find_contacts(const ContactContext &ctx, const PredA pred_a, const PredB pred_b, const search::SearchOptions opts, Tester &&tester) {

  using namespace search;
  using RecA = raw_predicate_arg_t<PredA>;
  using RecB = raw_predicate_arg_t<PredB>;

  static_assert(is_predicate_on<PredA, RecA>::value, "1st param must be a unary predicate returning a convertible to bool");
  static_assert(is_predicate_on<PredB, RecB>::value, "2nd param must be a unary predicate returning a convertible to bool");

  const auto &recs_a = ctx.ts.topo.records<RecA>();
  const auto &recs_b = ctx.ts.topo.records<RecB>();

  auto nr = neighbour_search<false>(search::CoordProvider<RecA, RecB>{ctx.ts}, recs_a, pred_a, recs_b, pred_b, opts, std::forward<Tester>(tester), ctx);

  return make_contacts<RecA, RecB>(ctx.ts.topo, nr);
}

} // namespace lahuta

#endif // LAHUTA_ENTITIES_FIND_CONTACTS_HPP
