#include "contacts/interactions.hpp"
#include "contacts/cationpi.hpp"
#include "contacts/halogen_bonds.hpp"
#include "contacts/hydrogen_bonds.hpp"
#include "contacts/hydrophobic.hpp"
#include "contacts/ionic.hpp"
#include "contacts/metalic.hpp"
#include "contacts/pistacking.hpp"

namespace lahuta {


Contacts Interactions::hbond() const { return find_hydrogen_bonds(luni_); }
Contacts Interactions::weak_hbond() const { return find_weak_hydrogen_bonds(luni_); }
Contacts Interactions::hydrophobic() const { return find_hydrophobic_bonds(luni_); }
Contacts Interactions::halogen() const { return find_halogen_bonds(luni_); }
Contacts Interactions::ionic() const { return find_ionic(luni_); }
Contacts Interactions::metalic() const { return find_metalic(luni_); }
Contacts Interactions::cationpi() const { return find_cationpi(luni_); }
Contacts Interactions::pistacking() const { return find_pistacking(luni_); }

} // namespace lahuta
