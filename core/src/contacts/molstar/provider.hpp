#ifndef LAHUTA_CONTACTS_MOLSTAR_PROVIDER_HPP
#define LAHUTA_CONTACTS_MOLSTAR_PROVIDER_HPP

#include "contacts/molstar/contacts.hpp"
#include "contacts/spec.hpp"
#include "entities/interaction_types.hpp"

namespace ms = lahuta::molstar;
// clang-format off
namespace lahuta {

struct MolStarContactProvider {
  ContactSpec<AtomRec,  AtomRec,  ms::HBondParams>       hbond       { InteractionType::HydrogenBond,      ms::make_hbond_recipe() };
  ContactSpec<AtomRec,  AtomRec,  ms::HBondParams>       whbond      { InteractionType::WeakHydrogenBond,  ms::make_weak_hbond_recipe() };
  ContactSpec<AtomRec,  AtomRec,  ms::HydrophobicParams> hydrophobic { InteractionType::Hydrophobic,       ms::make_hydrophobic_recipe() };
  ContactSpec<AtomRec,  AtomRec,  ms::HalogenParams>     halogen     { InteractionType::Halogen,           ms::make_halogen_recipe() };
  ContactSpec<AtomRec,  AtomRec,  ms::MetalicParams>     metalic     { InteractionType::MetalCoordination, ms::make_metalic_recipe() };
  ContactSpec<RingRec,  RingRec,  ms::PiStackingParams>  pistacking  { InteractionType::PiStacking,        ms::make_pistacking_recipe() };
  ContactSpec<GroupRec, RingRec,  ms::CationPiParams>    cationpi    { InteractionType::CationPi,          ms::make_cationpi_recipe() };
  ContactSpec<GroupRec, GroupRec, ms::IonicParams>       ionic       { InteractionType::Ionic,             ms::make_ionic_recipe() };

  constexpr auto specs() const noexcept {
    return std::tie(
      ionic, hbond, whbond, hydrophobic, halogen,
      pistacking, cationpi, metalic
    );
  }
};

} // namespace lahuta

#endif // LAHUTA_CONTACTS_MOLSTAR_PROVIDER_HPP
