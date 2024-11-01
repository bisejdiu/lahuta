#ifndef LAHUTA_HALOGEN_BONDS_HPP
#define LAHUTA_HALOGEN_BONDS_HPP

#include "atom_types.hpp"
#include "nn.hpp"

namespace lahuta {

class Luni;

inline struct HalogenParams {
  constexpr static double distance_max = 4.0;
  constexpr static double angle_max = M_PI / 6.0;                    // 30 degrees
  constexpr static double optimal_angle = M_PI;                      // 180 degrees
  constexpr static double optimal_acceptor_angle = M_PI * 2.0 / 3.0; // 120 degrees
} halogen_params;

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

Contacts find_halogen_bonds(const Luni &luni, HalogenParams opts = halogen_params);

} // namespace lahuta

#endif // LAHUTA_HALOGEN_BONDS_HPP
