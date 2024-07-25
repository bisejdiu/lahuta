#include <Geometry/point.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/RWMol.h>
#include <gemmi/mmread_gz.hpp>

using namespace gemmi;

RDKit::RWMol gemmiStructureToRDKit(Structure st, RDKit::Conformer &conf,
                                   bool ign_h = true);

RDKit::RWMol rdMolFromRDKitMol(RDKit::RWMol &mol,
                               std::vector<int> &atomIndices);
