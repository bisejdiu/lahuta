#ifndef LAHUTA_BONDS_PERCEPTION_SUBGRAPH_HPP
#define LAHUTA_BONDS_PERCEPTION_SUBGRAPH_HPP

#include <rdkit/GraphMol/RWMol.h>

#include "span.hpp"

namespace lahuta::bonds::subgraph {

// Indices must be strictly ascending and valid for source
RDKit::RWMol build_rdkit_submol(const RDKit::RWMol &source, span<const int> indices, bool include_bonds);

} // namespace lahuta::bonds::subgraph

#endif // LAHUTA_BONDS_PERCEPTION_SUBGRAPH_HPP
