#include "contacts/interactions.hpp"
#include "contacts/cationpi.hpp"
#include "contacts/hydrophobic.hpp"
#include "contacts/halogen_bonds.hpp"
#include "contacts/hydrogen_bonds.hpp"
#include "contacts/ionic.hpp"
#include "contacts/metalic.hpp"
#include "contacts/pistacking.hpp"

namespace lahuta {

ContactSet Interactions::hbond()        const {  return find_hydrogen_bonds(topology_); }
ContactSet Interactions::weak_hbond()   const {  return find_weak_hydrogen_bonds(topology_); }
ContactSet Interactions::hydrophobic()  const {  return find_hydrophobic_bonds(topology_); }
ContactSet Interactions::halogen()      const {  return find_halogen_bonds(topology_); }
ContactSet Interactions::ionic()        const {  return find_ionic(topology_); }
ContactSet Interactions::metalic()      const {  return find_metalic(topology_); }
ContactSet Interactions::cationpi()     const {  return find_cationpi(topology_); }
ContactSet Interactions::pistacking()   const {  return find_pistacking(topology_); }

} // namespace lahuta
