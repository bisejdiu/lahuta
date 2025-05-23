#ifndef LAHUTA_HALOGEN_BONDS_HPP
#define LAHUTA_HALOGEN_BONDS_HPP

#include "typing/types.hpp"
#include "entities/contact.hpp"
#include <GraphMol/Atom.h>
#include <GraphMol/RWMol.h>

namespace lahuta {
class Topology;

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

ContactSet find_halogen_bonds(const Topology &topology, const HalogenParams &params = HalogenParams{});

} // namespace lahuta

#endif // LAHUTA_HALOGEN_BONDS_HPP
