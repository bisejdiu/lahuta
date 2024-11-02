#include "contacts/interactions.hpp"
#include "contacts/cationpi.hpp"
#include "contacts/halogen_bonds.hpp"
#include "contacts/hydrogen_bonds.hpp"
#include "contacts/hydrophobic.hpp"
#include "contacts/ionic.hpp"
#include "contacts/metalic.hpp"
#include "contacts/pistacking.hpp"

namespace lahuta {

// FIX: these should functions should optionally take a custom function-specific options struct
// FIX: names: hydrophobic bonds are not really bonds, they are interactions or contacts

Contacts Interactions::find_hbond_interactions() { return find_hydrogen_bonds(luni_); }
Contacts Interactions::find_weak_hbond_interactions() { return find_weak_hydrogen_bonds(luni_); }
Contacts Interactions::find_hydrophobic_interactions() { return find_hydrophobic_bonds(luni_); }
Contacts Interactions::find_halogen_interactions() { return find_halogen_bonds(luni_); }
Contacts Interactions::find_ionic_interactions() { return find_ionic(luni_); }
Contacts Interactions::find_metalic_interactions() { return find_metalic(luni_); }
Contacts Interactions::find_cationpi_interactions() { return find_cationpi(luni_); }
Contacts Interactions::find_pistacking_interactions() { return find_pistacking(luni_); }

} // namespace lahuta
