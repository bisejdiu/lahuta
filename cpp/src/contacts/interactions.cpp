#include "contacts/interactions.hpp"
#include "contacts/distances.hpp"
#include "lahuta.hpp"

namespace lahuta {

Contacts Interactions::find_ionic_interactions() {
  PairFeatures finder(group_features_);
  finder.process_features(AtomType::POS_IONISABLE, AtomType::NEG_IONISABLE);

  Contacts interaction_contact = finder.find_interactions(InteractionType::Ionic, opts_.distance_cutoff);
  interaction_contact.set_luni(luni_);
  /*interaction_contact.print_interactions();*/
  return interaction_contact;
}

Contacts Interactions::find_hbond_interactions() {
  Contacts container(luni_);
  auto mol = luni_->get_molecule();
  const auto &atom_entities = luni_->get_atom_entities();
  auto atom_neighbors = luni_->find_neighbors2(6.0, 10);
  container.add_many(atom_neighbors, atom_entities);
  return container;
}

} // namespace lahuta
