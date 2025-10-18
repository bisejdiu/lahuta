#ifndef LAHUTA_CONTACTS_ARPEGGIO_CONTACTS_HPP
#define LAHUTA_CONTACTS_ARPEGGIO_CONTACTS_HPP

#include "contacts/recipe.hpp"
#include "entities/records.hpp"
#include "params.hpp"

namespace lahuta::arpeggio {

ContactRecipe<AtomRec,  AtomRec, HbondParams>           make_hbond_recipe();
ContactRecipe<AtomRec,  AtomRec, WeakHbondParams>       make_weak_hbond_recipe();
ContactRecipe<AtomRec,  AtomRec, PolarHbondParams>      make_polar_hbond_recipe();
ContactRecipe<AtomRec,  AtomRec, WeakPolarHbondParams>  make_weak_polar_hbond_recipe();
ContactRecipe<AtomRec,  AtomRec, HydrophobicParams>     make_hydrophobic_recipe();
ContactRecipe<AtomRec,  AtomRec, IonicParams>           make_ionic_recipe();
ContactRecipe<AtomRec,  AtomRec, AromaticParams>        make_aromatic_recipe();
ContactRecipe<AtomRec,  AtomRec, CarbonylParams>        make_carbonyl_recipe();
ContactRecipe<AtomRec,  AtomRec, VanDerWaalsParams>     make_vdw_recipe();
ContactRecipe<AtomRec,  RingRec, DonorPiParams>         make_donor_pi_recipe();
ContactRecipe<AtomRec,  RingRec, SulphurPiParams>       make_sulphur_pi_recipe();
ContactRecipe<AtomRec,  RingRec, CarbonPiParams>        make_carbon_pi_recipe();
ContactRecipe<AtomRec,  RingRec, CationPiParams>        make_cation_pi_recipe();

} // namespace lahuta::arpeggio

#endif // LAHUTA_CONTACTS_ARPEGGIO_CONTACTS_HPP
