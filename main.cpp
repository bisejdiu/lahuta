#include "Geometry/point.h"
#include "GraphMol/Atom.h"
#include "GraphMol/Bond.h"
#include "GraphMol/DetermineBonds/DetermineBonds.h"
#include "GraphMol/MonomerInfo.h"
#include "GraphMol/PeriodicTable.h"
#include "GraphMol/RWMol.h"
#include <iostream>

#include <gemmi/cif.hpp> // for Block
#include <gemmi/math.hpp>
// #include <gemmi/mmread.hpp> // for read_structure
// #include <gemmi/model.hpp> // for Model
// #include <gemmi/assembly.hpp>
#include <gemmi/mmread_gz.hpp> // for read_structure_gz
// #include <gemmi/neighbor.hpp>

#include "bonds.hpp"
#include "conv.hpp"
#include "nsgrid.hpp"

#include <chrono>

using namespace gemmi;

int main(int argc, char const *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <cif file>" << std::endl;
    return 1;
  }
  std::string file_name = argv[1];
  auto start = std::chrono::high_resolution_clock::now();
  Structure struc = read_structure_gz(file_name);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Time to load structure: " << elapsed.count() * 1000 << " ms"
            << std::endl;


  // std::map<std::string, std::vector<std::string>> atoms;
  // gemmi::cif::Block& block = gemmi::cif::read(file_name).sole_block();
  // auto atoms = block.find_mmcif_category("_atom_site.");
  // time rdKit conversion

  start = std::chrono::high_resolution_clock::now();
  RDKit::Conformer *conf = new RDKit::Conformer();
  RDKit::RWMol mol = gemmiStructureToRDKit(struc, *conf, false);
  mol.addConformer(conf, true);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to convert to RDKit: " << elapsed.count() * 1000 << " ms"
            << std::endl;

  // print atom information
  // RDKit::Atom *atom = mol.getAtomWithIdx(0);
  // RDKit::AtomPDBResidueInfo *res_info =
  //     (RDKit::AtomPDBResidueInfo *)atom->getMonomerInfo();
  // std::cout << "Atom name: " << atom->getSymbol() << std::endl;
  // std::cout << "Residue name: " << res_info->getResidueName() << std::endl;
  // std::cout << "Is hetero atom: " << res_info->getIsHeteroAtom() <<
  // std::endl;

  // // RDKit::determineConnectivity(mol, false, 0, 1.3, true); // 1271 bonds
  // // RDKit::determineConnectivity(mol, false, 0, 1.3, false); // 1269 bonds
  // // print number of models in the structure
  // std::cout << "Number of models: " << struc.models.size() << std::endl;
  //
  Model &model = struc.first_model();
  // time to find bonds
  start = std::chrono::high_resolution_clock::now();
  findBondsDeconstructed(struc, model, mol, 5.5, 1.3);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to find bonds using gemmi: " << elapsed.count() * 1000 << " ms"
            << std::endl;
  std::cout << "Number of bonds found using gemmi: " << mol.getNumBonds() << std::endl;
  //
  //
  // // print conformer information
  // RDKit::Conformer& conf2 = mol.getConformer();
  // std::cout << "Number of conformers: " << mol.getNumConformers() <<
  // std::endl; std::cout << "Number of atoms in conformer: " <<
  // conf2.getNumAtoms() << std::endl;
  //
  // // time for rdkit valence calculation
  start = std::chrono::high_resolution_clock::now();
  mol.updatePropertyCache(false);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to update property cache: " << elapsed.count() * 1000 << " ms"
            << std::endl;

  // iterate over all atoms
  // for (auto atom : mol.atoms()) {
  //
  //   // unsigned int dv =
  //   RDKit::PeriodicTable::getTable()->getDefaultValence(atom->getAtomicNum());
  //   // atom->calcExplicitValence(false);
  //
  //   std::cout
  //     << "Atom name: " << atom->getSymbol() << " "
  //     << atom->getExplicitValence() << std::endl;
  // }

  double cutoff = 5.0;

  start = std::chrono::high_resolution_clock::now();
  std::vector<RDGeom::Point3D> atom_coords = conf->getPositions();
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to extract all coords at once RDKit: " << elapsed.count() * 1000 << " ms"
            << std::endl;

  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to flatten coordinates: " << elapsed.count() * 1000 << " ms"
            << std::endl;

  // auto coords = flatten_coordinates(atom_coords);
  start = std::chrono::high_resolution_clock::now();
  FastNS grid(atom_coords, cutoff);
  auto results = grid.selfSearch();
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to find neighbors: " << elapsed.count() * 1000 << " ms"
            << std::endl;
  std::cout << "Number of neighbors: " << results.getNeighborPairsSize()
            << std::endl;

  start = std::chrono::high_resolution_clock::now();
  perceiveBonds(mol, results, 1.3);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to find bonds using NSGrid: " << elapsed.count() * 1000 << " ms"
            << std::endl;
  std::cout << "Number of bonds found using NSGrid: " << mol.getNumBonds() << std::endl;

  // start = std::chrono::high_resolution_clock::now();
  // grid.updateCutoff(5.5);
  // auto results2 = grid.selfSearch();
  // end = std::chrono::high_resolution_clock::now();
  // elapsed = end - start;
  // std::cout << "Time to find neighbors: " << elapsed.count() * 1000 << " ms"
  //           << std::endl;
  // std::cout << "Number of neighbors: " << results2.getNeighborPairsSize()
  //           << std::endl;
  //
  //
  // start = std::chrono::high_resolution_clock::now();
  // grid.updateCutoff(5.0);
  // auto results3 = grid.selfSearch();
  // end = std::chrono::high_resolution_clock::now();
  // elapsed = end - start;
  // std::cout << "Time to find neighbors: " << elapsed.count() * 1000 << " ms"
  //           << std::endl;
  // std::cout << "Number of neighbors: " << results3.getNeighborPairsSize()
  //           << std::endl;

  return 0;
}
