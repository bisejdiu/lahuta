#ifndef LAHUTA_MODEL_TOPOLOGY_HPP
#define LAHUTA_MODEL_TOPOLOGY_HPP

#include "GraphMol/Conformer.h"
#include "GraphMol/RWMol.h"
#include "models/parser.hpp"
#include <rdkit/GraphMol/MonomerInfo.h>

// clang-format off
using namespace RDKit;

namespace lahuta {

//
// Build Model Topology Under Fixed Assumptions:
// 1. Only standard protein amino acids
// 2. No H atoms, no missing atoms, no modified residues, no alternative locations
// 3. Atom ordering is consistent within the same residue
// 4. Only one chain (A)
// 5. Only one model
// 6. The last atom in the molecule is OXT
//

//
// Current Limitations:
// 1. Only CIF files are currently supported
// 2. No decompression is performed
//

//
// An important note about the design. By fixed assumptions, we mean that the state of the
// program is undefined if the input does not conform to these assumptions. This is a deliberate
// design choice, as it makes the code simpler and faster (no error checking, no branching, etc.).
// This also means that there are no, or very few, UB checks. The program from this point forward
// is in a "meaningful" state only as long as the assumptions are met.  - Besian, March 2025
//

enum class ModelTopologyMethod {
  None,
  Default,
  CSR
};

void build_model_topology_def(std::shared_ptr<RDKit::RWMol> &mol, RDKit::Conformer &conf, ModelParserResult &P);
void build_model_topology_csr(std::shared_ptr<RDKit::RWMol> &mol, RDKit::Conformer &conf, ModelParserResult &P);

inline void build_model_topology(std::shared_ptr<RDKit::RWMol> &mol, ModelParserResult &P, ModelTopologyMethod method = ModelTopologyMethod::Default) {
    auto conformer = std::make_unique<RDKit::Conformer>(mol->getNumAtoms());
    if (method == ModelTopologyMethod::CSR) {
       build_model_topology_csr(mol, *conformer.release(), P);
    } else if (method == ModelTopologyMethod::Default) {
        build_model_topology_def(mol, *conformer.release(), P);
    }
}

} // namespace lahuta

#endif // LAHUTA_MODEL_TOPOLOGY_HPP
