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
void read_and_build_model_topology(RDKit::RWMol &mol, RDKit::Conformer &conf, ModelParserResult &P);


// to help make atom info allocation more efficient
class AtomInfoContainer {
public:
  AtomInfoContainer(size_t count) {
    storage = static_cast<unsigned char*>(operator new[](count * sizeof(LeanAtomPDBResidueInfo)));
    capacity = count;
    constructed = 0;
  }

  ~AtomInfoContainer() {
    // Memory ownership will have been transferred to the Atom objects
  }

  /// Creates a new atom info object and transfers ownership
  std::unique_ptr<RDKit::AtomPDBResidueInfo>
  createAtomInfo(const char *atom_name, int serial, const char *res_name, int res_number) {
    if (constructed >= capacity) return nullptr; // Out of space

    // pointer to next available storage
    LeanAtomPDBResidueInfo *ptr = reinterpret_cast<LeanAtomPDBResidueInfo *>(&storage[constructed * sizeof(LeanAtomPDBResidueInfo)]);
    new (ptr) LeanAtomPDBResidueInfo(atom_name, serial, res_name, res_number); // in-place construction
    constructed++;

    return std::unique_ptr<RDKit::AtomPDBResidueInfo>(ptr);
  }

private:
  unsigned char *storage = nullptr;
  size_t capacity    = 0;
  size_t constructed = 0;

  AtomInfoContainer(const AtomInfoContainer &) = delete;
  AtomInfoContainer &operator=(const AtomInfoContainer &) = delete;
};


} // namespace lahuta

#endif // LAHUTA_MODEL_TOPOLOGY_HPP
