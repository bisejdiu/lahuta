#include "contacts/interactions.hpp"
#include "contacts/distances.hpp"
#include "contacts/halogen_bonds.hpp"
#include "contacts/hydrophobic.hpp"
#include "lahuta.hpp"

namespace lahuta {

// FIX: these should functions should optionally take a custom function-specific options struct
// FIX: atom neighbors are being computed per interaction type
// FIX: names: hydrophobic bonds are not really bonds, they are interactions or contacts
Contacts Interactions::find_ionic_interactions() {
  // FIX: Ionic contacts with 2ORG gives 229 contacts. which is likely wrong

  SimplePairFeatures finder(luni_->get_features());
  finder.process_features(AtomType::POS_IONISABLE, AtomType::NEG_IONISABLE);

  std::cout << "Distance cutoff: " << opts_.distance_cutoff << std::endl;
  Contacts interaction_contact = finder.find_interactions(InteractionType::Ionic, opts_.distance_cutoff);
  interaction_contact.set_luni(luni_);
  /*interaction_contact.print_interactions();*/
  return interaction_contact;
}

Contacts Interactions::find_hbond_interactions() {
  Contacts container(luni_);
  GeometryOptions _opts = GeometryOptions();
  find_hydrogen_bonds(*luni_, _opts, container);
  return container;
}

Contacts Interactions::find_weak_hbond_interactions() {
  Contacts container(luni_);
  GeometryOptions _opts = GeometryOptions();
  find_weak_hydrogen_bonds(*luni_, _opts, container);
  return container;
}

Contacts Interactions::find_hydrophobic_interactions() {
  Contacts container(luni_);
  GeometryOptions _opts = GeometryOptions();
  find_hydrophobic_bonds(*luni_, _opts, container);
  return container;
}

Contacts Interactions::find_halogen_interactions() {
  Contacts container(luni_);
  GeometryOptions _opts = GeometryOptions();
  find_halogen_bonds(*luni_, _opts, container);
  return container;
}

} // namespace lahuta
