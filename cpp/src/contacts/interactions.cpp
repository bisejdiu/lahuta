#include "contacts/interactions.hpp"
#include "contacts/cationpi.hpp"
#include "contacts/halogen_bonds.hpp"
#include "contacts/hydrogen_bonds.hpp"
#include "contacts/hydrophobic.hpp"
#include "contacts/ionic.hpp"
#include "contacts/metalic.hpp"
#include "contacts/pistacking.hpp"

namespace lahuta {


Contacts Interactions::hbond() { return find_hydrogen_bonds(luni_); }
Contacts Interactions::weak_hbond() { return find_weak_hydrogen_bonds(luni_); }
Contacts Interactions::hydrophobic() { return find_hydrophobic_bonds(luni_); }
Contacts Interactions::halogen() { return find_halogen_bonds(luni_); }
Contacts Interactions::ionic() { return find_ionic(luni_); }
Contacts Interactions::metalic() { return find_metalic(luni_); }
Contacts Interactions::cationpi() { return find_cationpi(luni_); }
Contacts Interactions::pistacking() { return find_pistacking(luni_); }

} // namespace lahuta
