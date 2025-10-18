#ifndef LAHUTA_MODEL_TOPOLOGY_HPP
#define LAHUTA_MODEL_TOPOLOGY_HPP

#include <rdkit/GraphMol/RWMol.h>

#include "models/parser.hpp"
#include "models/topology_impl.hpp"

// clang-format off
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
// - Only CIF files are supported
//

//
// An important note about the design. By fixed assumptions, we mean that the state of the
// program is undefined if the input does not conform to these assumptions. This is a deliberate
// design choice, as it makes the code simpler and faster (no error checking, no branching, etc.).
// This also means that there are no, or very few, UB checks. The program from this point forward
// is in a "meaningful" state only as long as the assumptions are met.  - Besian, March 2025
//
enum class ModelTopologyMethod { None, Default, CSR };

inline bool build_model_topology(std::shared_ptr<RDKit::RWMol> &mol, const ModelParserResult &P, ModelTopologyMethod method = ModelTopologyMethod::Default) {
    models::ModelTopology topology(P);

    models::ModelTopologyBuildingOptions options;
    if (method == ModelTopologyMethod::CSR) {
        options.graph_type = RDKit::GraphType::CSRMolGraph;
    } else {
        options.graph_type = RDKit::GraphType::MolGraph;
    }

    if (!topology.build(options)) {
        return false;
    }
    mol = topology.get_molecule();
    return true;
}

bool mock_build_model_topology(const ModelParserResult &P);

} // namespace lahuta

#endif // LAHUTA_MODEL_TOPOLOGY_HPP
