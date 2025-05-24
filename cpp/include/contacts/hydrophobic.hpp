#ifndef LAHUTA_HYDROPHOBIC_HPP
#define LAHUTA_HYDROPHOBIC_HPP

#include "entities/contact.hpp"
#include "typing/types.hpp"
#include <GraphMol/RWMol.h>

namespace lahuta {

class Topology;

struct HydrophobicParams {
  double distance_max = 4.0;
};

AtomType add_hydrophobic_atom(const RDKit::RWMol &mol, const RDKit::Atom &atom);

ContactSet find_hydrophobic_bonds(const Topology& topology, const HydrophobicParams& params = HydrophobicParams{});

} // namespace lahuta

#endif // LAHUTA_HYDROPHOBIC_HPP
