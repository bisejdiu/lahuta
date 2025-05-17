#ifndef LAHUTA_HYDROPHOBIC_HPP
#define LAHUTA_HYDROPHOBIC_HPP

#include "neighbors.hpp"
#include <GraphMol/RWMol.h>

namespace lahuta {

class Luni;

struct HydrophobicParams {
  double distance_max = 4.0;
};

AtomType add_hydrophobic_atom(const RDKit::RWMol &mol, const RDKit::Atom &atom);
Contacts find_hydrophobic_bonds(const Luni &luni, std::optional<HydrophobicParams> params = std::nullopt);

} // namespace lahuta

#endif // LAHUTA_HYDROPHOBIC_HPP
