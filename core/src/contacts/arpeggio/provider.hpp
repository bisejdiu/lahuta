#ifndef LAHUTA_CONTACTS_ARPEGGIO_PROVIDER_HPP
#define LAHUTA_CONTACTS_ARPEGGIO_PROVIDER_HPP

#include "contacts/arpeggio/contacts.hpp"
#include "contacts/spec.hpp"
#include "entities/interaction_types.hpp"

namespace ag = lahuta::arpeggio;

// clang-format off
namespace lahuta {

struct ArpeggioContactProvider {
  ContactSpec<AtomRec, AtomRec, ag::HbondParams>           hbond       { InteractionType::HydrogenBond,           ag::make_hbond_recipe() };
  ContactSpec<AtomRec, AtomRec, ag::WeakHbondParams>       whbond      { InteractionType::WeakHydrogenBond,       ag::make_weak_hbond_recipe() };
  ContactSpec<AtomRec, AtomRec, ag::PolarHbondParams>      phbond      { InteractionType::PolarHydrogenBond,      ag::make_polar_hbond_recipe() };
  ContactSpec<AtomRec, AtomRec, ag::WeakPolarHbondParams>  wphbond     { InteractionType::WeakPolarHydrogenBond,  ag::make_weak_polar_hbond_recipe() };
  ContactSpec<AtomRec, AtomRec, ag::HydrophobicParams>     hydrophobic { InteractionType::Hydrophobic,            ag::make_hydrophobic_recipe() };
  ContactSpec<AtomRec, AtomRec, ag::IonicParams>           ionic       { InteractionType::Ionic,                  ag::make_ionic_recipe() };
  ContactSpec<AtomRec, AtomRec, ag::AromaticParams>        aromatic    { InteractionType::Aromatic,               ag::make_aromatic_recipe() };
  ContactSpec<AtomRec, AtomRec, ag::CarbonylParams>        carbonyl    { InteractionType::Carbonyl,               ag::make_carbonyl_recipe() };
  ContactSpec<AtomRec, AtomRec, ag::VanDerWaalsParams>     vdw         { InteractionType::VanDerWaals,            ag::make_vdw_recipe() };
  ContactSpec<AtomRec, RingRec, ag::DonorPiParams>         donor_pi    { InteractionType::DonorPi,                ag::make_donor_pi_recipe() };
  ContactSpec<AtomRec, RingRec, ag::SulphurPiParams>       sulphur_pi  { InteractionType::SulphurPi,              ag::make_sulphur_pi_recipe() };
  ContactSpec<AtomRec, RingRec, ag::CarbonPiParams>        carbon_pi   { InteractionType::CarbonPi,               ag::make_carbon_pi_recipe() };
  ContactSpec<AtomRec, RingRec, ag::CationPiParams>        cation_pi   { InteractionType::CationPi,               ag::make_cation_pi_recipe() };

  constexpr auto specs() const noexcept {
    return std::tie(
      hbond, whbond, phbond, wphbond, hydrophobic,
      ionic, aromatic, carbonyl, vdw,
      donor_pi, sulphur_pi, carbon_pi, cation_pi
    );
  }
};

} // namespace lahuta

#endif // LAHUTA_CONTACTS_ARPEGGIO_PROVIDER_HPP

