/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   using Part = std::variant<const char*, std::string_view>;
 *   std::array<Part, 3> parts{Part{"besian"}, Part{"sejdiu"}, Part{"@gmail.com"}};
 *   std::string s;
 *   for (const auto& p : parts) {
 *     std::visit([&s](auto&& arg) { s += arg; }, p);
 *   }
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_CONTACTS_ENGINE_HPP
#define LAHUTA_CONTACTS_ENGINE_HPP

#include <optional>
#include <tuple>
#include <type_traits>
#include <utility>

#include "compute/topology_snapshot.hpp"
#include "entities/context.hpp"
#include "entities/find_contacts.hpp"
#include "recipe.hpp"
#include "spec.hpp"

namespace lahuta {
namespace C = lahuta::compute;

namespace {
template <typename T, typename = void>
struct has_specs : std::false_type {};

template <typename T>
struct has_specs<T, std::void_t<decltype(std::declval<const T>().specs())>> : std::true_type {};

template <typename>
struct is_contact_spec : std::false_type {};

template <typename RecA, typename RecB, typename Params>
struct is_contact_spec<ContactSpec<RecA, RecB, Params>> : std::true_type {};

template <typename>
struct all_contact_specs : std::false_type {};
template <typename... Ts>
struct all_contact_specs<std::tuple<Ts...>> : std::conjunction<is_contact_spec<std::decay_t<Ts>>...> {};
} // namespace

template <class Provider>
class InteractionEngine : public Provider {
  static_assert(has_specs<Provider>::value, "Provider must define a constexpr specs() method");
  using specs_t = decltype(std::declval<const Provider>().specs());
  static_assert(std::is_same_v<specs_t, std::remove_cv_t<specs_t>>, "specs() must return a std::tuple<T...>");
  static_assert(all_contact_specs<specs_t>::value,
                "All element of Provider::specs() must be ContactSpec<Rec,Rec,Params>");

public:
  using Provider::Provider;

  [[nodiscard]] ContactSet compute(const C::TopologySnapshot &ts) const {
    return compute(ts, /*no tag filter*/ std::nullopt);
  }

  [[nodiscard]] ContactSet compute(const C::TopologySnapshot &ts,
                                   std::optional<InteractionTypeSet> only) const {
    ContactSet out;
    auto run_spec = [&](auto &spec) {
      if (!spec.enabled) return;
      if (only && !only->contains(spec.tag)) return;

      Logger::get_logger()->debug("InteractionEngine: Computing {} contacts",
                                  interaction_type_to_string(spec.tag));

      auto &R   = spec.recipe;
      auto ctx  = ContactContext{ts, R.params};
      auto opts = search::make_search_opts(spec.tag.category);

      opts.distance_max = R.params.distance_max; // we need to override the default distance_max, which is not
                                                 // guaranteed to be set or correct

      if (R.mode == RecipeMode::Self) {
        out.insert(find_contacts(ctx, R.pred_a, opts, R.tester));
      } else {
        out.insert(find_contacts(ctx, R.pred_a, R.pred_b, opts, R.tester));
      }
      Logger::get_logger()->debug("InteractionEngine: Completed {} contacts computation",
                                  interaction_type_to_string(spec.tag));
    };

    std::apply([&](auto &...Ss) { (run_spec(Ss), ...); }, this->specs());
    return out;
  }

  void set_enabled(InteractionType tag, bool on = true) {
    auto flip = [&](auto &spec) {
      if (spec.tag == tag) spec.enabled = on;
    };
    std::apply([&](auto &...Ss) { (flip(Ss), ...); }, this->specs());
  }
};

} // namespace lahuta

#endif // LAHUTA_CONTACTS_ENGINE_HPP
