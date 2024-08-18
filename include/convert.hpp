#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>
#include <gemmi/mmread_gz.hpp>

using namespace gemmi;

void gemmiStructureToRDKit(RDKit::RWMol &mol, const Structure &st, RDKit::Conformer &conf,
                                   bool ign_h = true);

RDKit::RWMol rdMolFromRDKitMol(RDKit::RWMol &mol, std::vector<int> &atomIndices);
// RDKit::RWMol _rdMolFromRDKitMol(RDKit::RWMol &mol, std::vector<int> &atomIndices);
// RDKit::RWMol rdMolFromRDKitMol(RDKit::RWMol &mol, std::vector<std::optional<int>> &atomIndices);


RDKit::RWMol filter_atoms(RDKit::RWMol &mol, std::vector<int> &atomIndices);
RDKit::RWMol filter_atom_conf(RDKit::RWMol &mol, std::vector<int> &atomIndices);
RDKit::RWMol filter_with_atom_data(RDKit::RWMol &mol, std::vector<int> &atomIndices);

