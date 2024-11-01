#ifndef LAHUTA_METALS_HPP
#define LAHUTA_METALS_HPP

#include "GraphMol/RWMol.h"
#include "atom_types.hpp"

namespace lahuta {

// clang-format off
inline std::vector<int> IonicTypeMetals = {
    3, 11, 19, 37, 55,  // Li, Na, K, Rb, Cs
    12, 20, 38, 56,     // Mg, Ca, Sr, Ba
    13, 31, 49, 81,     // Al, Ga, In, Tl
    21, 50, 82, 83, 51, // Sc, Sn, Pb, Bi, Sb
    80                  // Hg
};
// clang-format on

bool is_halogen(int atomic_num);

bool is_protein_sidechain(const std::string &atomname);

bool is_protein_backbone(const std::string &atomname);

bool is_nucleic_backbone(const std::string &atomname);

bool is_transition_metal(int atomic_num);

AtomType add_metal(const RDKit::RWMol &mol, const RDKit::Atom &atom);

AtomType add_metal_binding(const RDKit::RWMol &mol, const RDKit::Atom &atom);

} // namespace lahuta

#endif // LAHUTA_METALS_HPP
