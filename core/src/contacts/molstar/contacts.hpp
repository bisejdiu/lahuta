/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s = "besian@gmail.com";
 *   s.insert(6, "sejdiu");
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_CONTACTS_MOLSTAR_CONTACTS_HPP
#define LAHUTA_CONTACTS_MOLSTAR_CONTACTS_HPP

#include "contacts/recipe.hpp"
#include "entities/records.hpp"
#include "params.hpp"

// clang-format off
namespace lahuta::molstar {

ContactRecipe<AtomRec,  AtomRec, HBondParams>       make_hbond_recipe();
ContactRecipe<AtomRec,  AtomRec, HBondParams>       make_weak_hbond_recipe();
ContactRecipe<AtomRec,  AtomRec, HydrophobicParams> make_hydrophobic_recipe();
ContactRecipe<AtomRec,  AtomRec, HalogenParams>     make_halogen_recipe();
ContactRecipe<AtomRec,  AtomRec, MetalicParams>     make_metalic_recipe();
ContactRecipe<GroupRec, RingRec, CationPiParams>    make_cationpi_recipe();
ContactRecipe<RingRec,  RingRec, PiStackingParams>  make_pistacking_recipe();
ContactRecipe<GroupRec, GroupRec,IonicParams>       make_ionic_recipe();

} // namespace lahuta::molstar

#endif // LAHUTA_CONTACTS_MOLSTAR_CONTACTS_HPP
