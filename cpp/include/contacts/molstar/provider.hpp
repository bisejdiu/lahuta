#ifndef LAHUTA_CONTACTS_MOLSTAR_PROVIDER_HPP
#define LAHUTA_CONTACTS_MOLSTAR_PROVIDER_HPP

#include "contacts/molstar/contacts.hpp"
#include "contacts/spec.hpp"

// clang-format off
namespace lahuta {

struct MolStarContactProvider {
  ContactSpec<AtomRec,  AtomRec,  HBondParameters>   hbond       { InteractionType::HydrogenBond,      make_hbond_recipe() };
  ContactSpec<AtomRec,  AtomRec,  HBondParameters>   whbond      { InteractionType::WeakHydrogenBond,  make_weak_hbond_recipe() };
  ContactSpec<AtomRec,  AtomRec,  HydrophobicParams> hydrophobic { InteractionType::Hydrophobic,       make_hydrophobic_recipe() };
  ContactSpec<AtomRec,  AtomRec,  HalogenParams>     halogen     { InteractionType::Halogen,           make_halogen_recipe() };
  ContactSpec<AtomRec,  AtomRec,  MetalicParams>     metalic     { InteractionType::MetalCoordination, make_metalic_recipe() };
  ContactSpec<RingRec,  RingRec,  PiStackingParams>  pistacking  { InteractionType::PiStacking,        make_pistacking_recipe() };
  ContactSpec<GroupRec, RingRec,  CationPiParams>    cationpi    { InteractionType::CationPi,          make_cationpi_recipe() };
  ContactSpec<GroupRec, GroupRec, IonicParams>       ionic       { InteractionType::Ionic,             make_ionic_recipe() };

  constexpr auto specs() const noexcept {
    return std::tie(
      ionic, hbond, whbond, hydrophobic, halogen,
      pistacking, cationpi, metalic
    );
  }
};

} // namespace lahuta

#endif // LAHUTA_CONTACTS_MOLSTAR_PROVIDER_HPP
