#ifndef LAHUTA_ENTITIES_FIND_CONTACTS_HPP
#define LAHUTA_ENTITIES_FIND_CONTACTS_HPP

#include "contact.hpp"
#include "context.hpp"
#include "entity_id.hpp"
#include "pipeline/function_traits.hpp"
#include "records.hpp"
#include "search/hit_buffer.hpp"
#include "search/provider.hpp"
#include "search/search.hpp"
#include "topology.hpp"
#include <type_traits>
#include <vector>

// clang-format off
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
template<typename F>
using raw_predicate_arg_t = std::decay_t<util::function_arg_t<std::decay_t<F>, 0>>;

template<typename F, typename Rec>
using is_predicate_on = std::integral_constant<bool,
  std::is_invocable<F, const Rec&>::value &&
  std::is_convertible<std::invoke_result_t<F, const Rec&>, bool>::value>;

template<
  typename Pred,
  typename Tester,
  typename Rec = raw_predicate_arg_t<Pred>,
  std::enable_if_t<is_predicate_on<Pred, Rec>::value &&
                   std::is_invocable_r_v<InteractionType, Tester, uint32_t, uint32_t, float, const ContactContext&>, int> = 0>
ContactSet find_contacts(const ContactContext& ctx, const Pred pred, const search::SearchOptions opts, Tester&& tester) {
  using namespace search;

  const auto &recs = ctx.topology.records<Rec>();
  auto result = neighbour_search<true>(CoordProvider<Rec, Rec>{ctx.topology}, recs, pred, recs, pred, opts, std::forward<Tester>(tester), ctx);

  return make_contacts<Rec, Rec>(ctx.topology, result);
}

template<
  typename PredA, typename PredB, typename Tester,
  std::enable_if_t<std::is_invocable_r_v<InteractionType, Tester, uint32_t, uint32_t, float, const ContactContext&>, int> = 0>
ContactSet find_contacts(const ContactContext& ctx, const PredA pred_a, const PredB pred_b, const search::SearchOptions opts, Tester &&tester) {

  using namespace search;
  using RecA = raw_predicate_arg_t<PredA>;
  using RecB = raw_predicate_arg_t<PredB>;

  static_assert(is_predicate_on<PredA, RecA>::value, "1st param must be a unary predicate returning bool‑convertible");
  static_assert(is_predicate_on<PredB, RecB>::value, "2nd param must be a unary predicate returning bool‑convertible");

  const auto &recs_a = ctx.topology.records<RecA>();
  const auto &recs_b = ctx.topology.records<RecB>();

  auto nr = neighbour_search<false>(CoordProvider<RecA, RecB>{ctx.topology}, recs_a, pred_a, recs_b, pred_b, opts, std::forward<Tester>(tester), ctx);

  return make_contacts<RecA, RecB>(ctx.topology, nr);
}

} // namespace lahuta

#endif // LAHUTA_ENTITIES_FIND_CONTACTS_HPP
