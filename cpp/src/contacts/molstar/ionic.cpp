#include "contacts/molstar/contacts.hpp"

// clang-format off
namespace lahuta {

ContactRecipe<GroupRec, GroupRec, IonicParams> make_ionic_recipe() {
  return {
    IonicParams{},
    +[](GroupRec const& r){ return (r.a_type & AtomType::PositiveCharge) == AtomType::PositiveCharge; },
    +[](GroupRec const& r){ return (r.a_type & AtomType::NegativeCharge) == AtomType::NegativeCharge; },
    +[](u32 a, u32 b, float d, ContactContext const&){
      if (d < 2.0f || a == b) return InteractionType::None;
      return InteractionType::Ionic;
    }
  };
}

} // namespace lahuta
