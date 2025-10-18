#ifndef LAHUTA_HYDROPHOBIC_HPP
#define LAHUTA_HYDROPHOBIC_HPP

#include <GraphMol/RWMol.h>

#include "typing/types.hpp"

namespace lahuta {

AtomType add_hydrophobic_atom(const RDKit::RWMol &mol, const RDKit::Atom &atom);

} // namespace lahuta

#endif // LAHUTA_HYDROPHOBIC_HPP
