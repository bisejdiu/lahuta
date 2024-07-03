#include "Geometry/point.h"
#include "GraphMol/Atom.h"
#include "GraphMol/Bond.h"
#include "GraphMol/DetermineBonds/DetermineBonds.h"
#include "GraphMol/MonomerInfo.h"
#include "GraphMol/PeriodicTable.h"
#include "GraphMol/RDKitBase.h"
#include "GraphMol/RWMol.h"
#include <iostream>

#include <gemmi/cif.hpp> // for Block
#include <gemmi/math.hpp>
#include <gemmi/mmread_gz.hpp> // for read_structure_gz
#include <gemmi/model.hpp>     // for Model

#include "bonds.hpp"
#include "conv.hpp"
#include "nsgrid.hpp"

#include <chrono>
#include <gemmi/metadata.hpp>

#include "bond_utils.hpp"

using namespace gemmi;


int main(int argc, char const *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <cif file>" << std::endl;
    return 1;
  }
  std::string file_name = argv[1];
  auto start = std::chrono::high_resolution_clock::now();
  Structure st = read_structure_gz(file_name);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Time to load structure: " << elapsed.count() * 1000 << " ms"
            << std::endl;

  start = std::chrono::high_resolution_clock::now();
  RDKit::Conformer *conf = new RDKit::Conformer();
  RDKit::RWMol mol = gemmiStructureToRDKit(st, *conf, false);
  mol.addConformer(conf, true);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to convert to RDKit: " << elapsed.count() * 1000 << " ms"
            << std::endl;

  // update property cache
  // mol.updatePropertyCache(false);

  double cutoff = 4.5; // 01;
  std::vector<RDGeom::Point3D> atom_coords = conf->getPositions();
  FastNS grid(atom_coords, cutoff);
  auto results = grid.selfSearch();
  std::cout << "Number of neighbors: " << results.getNeighborPairsSize()
            << std::endl;

  // perceives bonds
  start = std::chrono::high_resolution_clock::now();
  perceiveBonds(mol, results, 0.45);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to perceive bonds: " << elapsed.count() * 1000 << " ms"
            << std::endl;
  std::cout << "Number of bonds found using NSGrid: " << mol.getNumBonds()
            << std::endl;

  // mol.updatePropertyCache(false);
  // compare atom.getExplicitValence() with GetExplicitValenceFromAtomBonds
  // NOTE: the valences do not differ, but RDKit atoms & bonds need updating,
  // before getting the correct valence for (const auto &atom: mol.atoms()) {
  //   auto valence = GetExplicitValenceFromAtomBonds(mol, atom);
  //   if (atom->getExplicitValence() != valence) {
  //     std::cout << "V Diff: Atom: " << atom->getIdx() << " Valence: " <<
  //     atom->getExplicitValence() << " " << valence << std::endl;
  //   }
  // }

  // Clean up incorrect bonds
  start = std::chrono::high_resolution_clock::now();
  CleanUpMolecule(mol, *conf);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to clean up molecule: " << elapsed.count() * 1000 << " ms"
            << std::endl;

  for (gemmi::Connection &conn : st.connections) {
    // FIX: Need to iterate over all models
    Atom *a1 = st.first_model().find_cra(conn.partner1).atom;
    Atom *a2 = st.first_model().find_cra(conn.partner2).atom;

    if (mol.getBondBetweenAtoms(a1->serial - 1, a2->serial - 1) == nullptr) {
      // std::cout << "Adding bond: " << a1->serial << " - " << a2->serial <<
      // std::endl;
      mol.addBond((unsigned int)a1->serial - 1, (unsigned int)a2->serial - 1,
                  RDKit::Bond::BondType::SINGLE);
      continue;
    }
    // std::cout << "Bond already exists: " << a1->serial << " - " << a2->serial
    // << std::endl;
  }

  // Clean up incorrect bonds
  // start = std::chrono::high_resolution_clock::now();
  // CleanUpMolecule(mol, *conf);
  // end = std::chrono::high_resolution_clock::now();
  // elapsed = end - start;
  // std::cout << "Time to clean up molecule: " << elapsed.count() * 1000 << "
  // ms"
  //           << std::endl;

  // update property cache AGAIN
  mol.updatePropertyCache(false);

  // (Should we) Clean up incorrect bonds again?
  // CleanUpMolecule(mol, *conf);

  std::cout << "Number of bonds w/ _struct_conns " << mol.getNumBonds()
            << std::endl;

  return 0;
}
