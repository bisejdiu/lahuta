#ifndef LAHUTA_CONTACTS_RECIPE_HPP
#define LAHUTA_CONTACTS_RECIPE_HPP

#include <cstdint>

#include "entities/context.hpp"
#include "entities/entity_id.hpp"
#include "entities/interaction_types.hpp"
#include "entities/records.hpp"

// clang-format off
namespace lahuta {

// enum class RecipeState : uint8_t { Disabled, Enabled };
enum class RecipeMode  : uint8_t { Self, Cross };

template<typename RecA, typename RecB, typename Params>
struct ContactRecipe {
  using rec_a_t = RecA;
  using rec_b_t = RecB;

  static constexpr Kind kind_a = KindOf<RecA>::value;
  static constexpr Kind kind_b = KindOf<RecB>::value;

  Params params {};
  RecipeMode mode = RecipeMode::Cross;

  bool (*pred_a)(const RecA&) = nullptr;
  bool (*pred_b)(const RecB&) = nullptr;
  InteractionType (*tester)(std::uint32_t,std::uint32_t, float,const ContactContext&) = nullptr;

  template<typename PredA, typename PredB, typename Tester>
  constexpr ContactRecipe(Params p, PredA pa, PredB pb, Tester ts ) noexcept
    : params{p}, pred_a{+pa}, pred_b{+pb}, tester{+ts} {
      static_assert(std::is_convertible_v<PredA, bool (*)(const RecA&)>, "PredA must be captureless");
      static_assert(std::is_convertible_v<PredB, bool (*)(const RecB&)>, "PredB must be captureless");
      static_assert(std::is_convertible_v<Tester, InteractionType (*)(std::uint32_t,std::uint32_t, float,const ContactContext&)>, "Tester must be captureless");
    }

  template<typename Pred, typename Tester>
  constexpr ContactRecipe(Params p, Pred pred, Tester ts ) noexcept
    : params{p}, pred_a{+pred}, pred_b{+pred}, tester{+ts}, mode{RecipeMode::Self} {
        static_assert(std::is_convertible_v<Pred, bool (*)(const RecA&)>, "Pred must be captureless");
        static_assert(std::is_convertible_v<Tester, InteractionType (*)(std::uint32_t,std::uint32_t, float,const ContactContext&)>, "Tester must be captureless");
        static_assert(std::is_same_v<RecA,RecB>, "One-predicate constructor is only valid when RecA == RecB");
    }
};

template<typename RecA, typename RecB, typename Params>
ContactRecipe(Params, bool (*)(const RecA&), bool (*)(const RecB&), InteractionType (*)(u32,u32,float,const ContactContext&)) -> ContactRecipe<RecA,RecB,Params>;

template<typename Rec, typename Params>
ContactRecipe(Params, bool (*)(const Rec&),                         InteractionType (*)(u32,u32,float,const ContactContext&)) -> ContactRecipe<Rec,Rec,Params>;

} // namespace lahuta

#endif // LAHUTA_CONTACTS_RECIPE_HPP
