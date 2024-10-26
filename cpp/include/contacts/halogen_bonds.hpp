#ifndef LAHUTA_HALOGEN_BONDS_HPP
#define LAHUTA_HALOGEN_BONDS_HPP

#include "atom_types.hpp"
#include "contacts/hydrogen_bonds.hpp"
#include "nn.hpp"

namespace lahuta {

class Luni;

inline std::unordered_set<int> hal_bond_elements = {17, 35, 53};
inline std::unordered_set<int> X = {7, 8, 16};
inline std::unordered_set<int> Y = {6, 7, 15, 16};

/**
 * Halogen bond donors (X-C, with X one of Cl, Br, I or At) not F!
 */
AtomType add_halogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom);

/**
 * Halogen bond acceptors (Y-{O|N|S}, with Y=C,P,N,S)
 */
AtomType add_halogen_acceptor(const RDKit::RWMol &mol, const RDKit::Atom &atom);

void find_halogen_bonds(Luni &luni, const GeometryOptions &opts, Contacts &container);

struct HalogenBondsParams {
  double distanceMax = 4.0; // Maximum distance for halogen bonds
  double angleMax = 30.0;   // Maximum angle for halogen bonds
};

struct HalogenBondsOptions {
  double angleMax; // in radians
};

} // namespace lahuta

#endif // LAHUTA_HALOGEN_BONDS_HPP
