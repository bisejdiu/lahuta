#ifndef LAHUTA_CONTACTS_GETCONTACTS_CONTACTS_HPP
#define LAHUTA_CONTACTS_GETCONTACTS_CONTACTS_HPP

#include "contacts/getcontacts/params.hpp"
#include "contacts/recipe.hpp"
#include "entities/records.hpp"

// clang-format off
namespace lahuta::getcontacts {

ContactRecipe<AtomRec, AtomRec, HydrogenBondParams> make_hbond_recipe();
ContactRecipe<AtomRec, AtomRec, HydrophobicParams>  make_hydrophobic_recipe();
ContactRecipe<AtomRec, AtomRec, SaltBridgeParams>   make_salt_bridge_recipe();
ContactRecipe<AtomRec, AtomRec, VanDerWaalsParams>  make_vdw_recipe();
ContactRecipe<AtomRec, RingRec, PiCationParams>     make_pi_cation_recipe();
ContactRecipe<RingRec, RingRec, PiStackingParams>   make_pi_stacking_recipe();
ContactRecipe<RingRec, RingRec, TStackingParams>    make_t_stacking_recipe();

} // namespace lahuta::getcontacts

#endif // LAHUTA_CONTACTS_GETCONTACTS_CONTACTS_HPP
