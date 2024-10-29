#include "contacts/interactions.hpp"
#include "contacts/distances.hpp"
#include "contacts/halogen_bonds.hpp"
#include "contacts/hydrophobic.hpp"
#include "contacts/ionic.hpp"
#include "contacts/metalic.hpp"
#include "lahuta.hpp"

namespace lahuta {

// FIX: Luni is sometimes passed as a pointer, sometimes as a reference (should be consistent)
// FIX: these should functions should optionally take a custom function-specific options struct
// FIX: atom neighbors are being computed per interaction type
// FIX: names: hydrophobic bonds are not really bonds, they are interactions or contacts

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

Contacts Interactions::find_ionic_interactions() {
  // FIX: Since we have to pass Luni to the finder, we can just have it return the contacts
  // directly, making it unnecessary for us to instantiate a Contacts object (which also takes a
  // Luni object)
  Contacts contacts(luni_);
  GeometryOptions _opts_;
  find_ionic(*luni_, _opts_, contacts);
  return contacts;
}

Contacts Interactions::find_metalic_interactions() {
  Contacts contacts(luni_);
  GeometryOptions _opts_;
  find_metalic(luni_, _opts_, contacts);
  return contacts;
}

} // namespace lahuta
