#ifndef LAHUTA_HALOGEN_BONDS_HPP
#define LAHUTA_HALOGEN_BONDS_HPP

#include "neighbors.hpp"

namespace lahuta {

class Luni;

struct HalogenParams {
  double distance_max = 4.0;
  double angle_max = M_PI / 6.0;                    // 30 degrees
  double optimal_angle = M_PI;                      // 180 degrees
  double optimal_acceptor_angle = M_PI * 2.0 / 3.0; // 120 degrees
};

// Halogen bond donors (X-C, with X one of Cl, Br, I or At) not F!
AtomType add_halogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom);

// Halogen bond acceptors (Y-{O|N|S}, with Y=C,P,N,S)
AtomType add_halogen_acceptor(const RDKit::RWMol &mol, const RDKit::Atom &atom);

Contacts find_halogen_bonds(const Luni &luni, std::optional<HalogenParams> params = std::nullopt);

} // namespace lahuta

#endif // LAHUTA_HALOGEN_BONDS_HPP
