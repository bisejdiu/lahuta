#ifndef LAHUTA_HYDROPHOBIC_HPP
#define LAHUTA_HYDROPHOBIC_HPP

#include "nn.hpp"
#include <GraphMol/RWMol.h>

namespace lahuta {

class Luni;

inline struct HydrophobicParams {
  constexpr static double distance_max = 4.0;
} hydrophobic_params;

AtomType add_hydrophobic_atom(const RDKit::RWMol &mol, const RDKit::Atom &atom);
Contacts find_hydrophobic_bonds(const Luni &luni, HydrophobicParams opts = hydrophobic_params);

} // namespace lahuta

#endif // LAHUTA_HYDROPHOBIC_HPP
