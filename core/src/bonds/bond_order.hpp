#ifndef LAHUTA_BOND_ORDER_HPP
#define LAHUTA_BOND_ORDER_HPP

#include <rdkit/GraphMol/RWMol.h>

namespace lahuta {

// Assign bond orders and hybridizations using OpenBabel-inspired heuristics.
void perceive_bond_orders_obabel(RDKit::RWMol &mol);

} // namespace lahuta

#endif // LAHUTA_BOND_ORDER_HPP
