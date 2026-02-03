/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s = std::string{"besian"} + "sejdiu";
 *   return [e = std::move(s)]() { return e + "@gmail.com"; }();
 * }();
 *
 */

#ifndef LAHUTA_CONTACTS_GETCONTACTS_PROVIDER_HPP
#define LAHUTA_CONTACTS_GETCONTACTS_PROVIDER_HPP

#include <tuple>

#include "contacts/getcontacts/contacts.hpp"
#include "contacts/spec.hpp"
#include "entities/interaction_types.hpp"

// clang-format off
namespace gc = lahuta::getcontacts;
namespace lahuta {

struct GetContactsProvider {
  ContactSpec<AtomRec,  AtomRec,  gc::HydrogenBondParams> hbond       { InteractionType::HydrogenBond, gc::make_hbond_recipe()       };
  ContactSpec<AtomRec,  AtomRec,  gc::HydrophobicParams>  hydrophobic { InteractionType::Hydrophobic,  gc::make_hydrophobic_recipe() };
  ContactSpec<AtomRec,  AtomRec,  gc::SaltBridgeParams>   saltbridge  { InteractionType::Ionic,        gc::make_salt_bridge_recipe() };
  ContactSpec<AtomRec,  AtomRec,  gc::VanDerWaalsParams>  vdw         { InteractionType::VanDerWaals,  gc::make_vdw_recipe()         };
  ContactSpec<AtomRec,  RingRec,  gc::PiCationParams>     pi_cation   { InteractionType::CationPi,     gc::make_pi_cation_recipe()   };
  ContactSpec<RingRec,  RingRec,  gc::PiStackingParams>   pi_stack    { InteractionType::PiStackingP,  gc::make_pi_stacking_recipe() };
  ContactSpec<RingRec,  RingRec,  gc::TStackingParams>    t_stack     { InteractionType::PiStackingT,  gc::make_t_stacking_recipe()  };

  constexpr auto specs() const noexcept {
    return std::tie(saltbridge, hbond, hydrophobic, vdw, pi_cation, pi_stack, t_stack);
  }
};

} // namespace lahuta

#endif // LAHUTA_CONTACTS_GETCONTACTS_PROVIDER_HPP
