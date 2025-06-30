#ifndef LAHUTA_CONTACTS_ENGINE_HPP
#define LAHUTA_CONTACTS_ENGINE_HPP

#include "entities/context.hpp"
#include "entities/find_contacts.hpp"
#include "recipe.hpp"
#include "spec.hpp"

// clang-format off
namespace lahuta {

namespace {
template<typename T, typename=void>
struct has_specs : std::false_type {};

template<typename T>
struct has_specs<T, std::void_t<decltype(std::declval<const T>().specs())>> : std::true_type {};

template<typename> struct is_contact_spec : std::false_type {};

template<typename RecA, typename RecB, typename Params>
struct is_contact_spec<ContactSpec<RecA,RecB,Params>> : std::true_type {};

template<typename> struct all_contact_specs : std::false_type {};
template<typename... Ts>
struct all_contact_specs<std::tuple<Ts...>> : std::conjunction<is_contact_spec<std::decay_t<Ts>>...> {};
} // namespace

template<class Provider>
class InteractionEngine : public Provider {
  static_assert(has_specs<Provider>::value, "Provider must define a constexpr specs() method");
  using specs_t = decltype(std::declval<const Provider>().specs());
  static_assert(std::is_same_v<specs_t, std::remove_cv_t<specs_t>>, "specs() must return a std::tuple<T...>");
  static_assert(all_contact_specs<specs_t>::value, "All element of Provider::specs() must be ContactSpec<Rec,Rec,Params>");

public:
  using Provider::Provider;

  [[nodiscard]] ContactSet compute(const Topology& topo) const {
    return compute(topo, /*no tag filter*/ std::nullopt);
  }

  [[nodiscard]] ContactSet compute(const Topology& topo, std::optional<InteractionType> only) const {

    ContactSet out;
    auto run_spec = [&](auto& spec){
      if (!spec.enabled) return;
      if (only && spec.tag != *only) return;

      auto& R   = spec.recipe;
      auto ctx  = ContactContext{topo, R.params};
      auto opts = search::make_search_opts(spec.tag.category);
      opts.distance_max = R.params.distance_max; // we need to override the default distance_max, which is not guaranteed to be set or correct

      if (R.mode == RecipeMode::Self) {
        out.insert(find_contacts(ctx, R.pred_a, opts, R.tester));
      } else {
        out.insert(find_contacts(ctx, R.pred_a, R.pred_b, opts, R.tester));
      }
    };

    std::apply([&](auto&... Ss){ (run_spec(Ss), ...); }, this->specs());
    return out;
  }

  void set_enabled(InteractionType tag, bool on=true) {
    auto flip = [&](auto& spec) {if (spec.tag == tag) spec.enabled = on;};
    std::apply([&](auto&... Ss){ (flip(Ss), ...); }, this->specs());
  }
};

} // namespace lahuta

#endif // LAHUTA_CONTACTS_ENGINE_HPP
