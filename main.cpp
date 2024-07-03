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

constexpr const char *hybridizationToString(
    RDKit::Atom::HybridizationType type) {
  switch (type) {
    case RDKit::Atom::HybridizationType::UNSPECIFIED:
      return "";
    case RDKit::Atom::HybridizationType::S:
      return "S";
    case RDKit::Atom::HybridizationType::SP:
      return "SP";
    case RDKit::Atom::HybridizationType::SP2:
      return "SP2";
    case RDKit::Atom::HybridizationType::SP3:
      return "SP3";
    case RDKit::Atom::HybridizationType::SP3D:
      return "SP3D";
    case RDKit::Atom::HybridizationType::SP2D:
      return "SP2D";
    case RDKit::Atom::HybridizationType::SP3D2:
      return "SP3D2";
    case RDKit::Atom::HybridizationType::OTHER:
      return "OTHER";
  }
  return "";
}

double AverageBondAngle(const RDKit::ROMol &mol, const RDKit::Atom *atom,
                        const RDKit::Conformer &conf) {
  double avgDegrees = 0.0;
  int n = 0;

  for (auto bondIt1 = mol.getAtomBonds(atom); bondIt1.first != bondIt1.second;
       ++bondIt1.first) {
    const RDKit::Bond *bond1 = mol[*bondIt1.first];
    const RDKit::Atom *neighbor1 = bond1->getOtherAtom(atom);

    for (auto bondIt2 = bondIt1; ++bondIt2.first != bondIt2.second;) {
      const RDKit::Bond *bond2 = mol[*bondIt2.first];
      const RDKit::Atom *neighbor2 = bond2->getOtherAtom(atom);

      // Get the coordinates of the atoms
      RDGeom::Point3D pos1 = conf.getAtomPos(atom->getIdx());
      RDGeom::Point3D pos2 = conf.getAtomPos(neighbor1->getIdx());
      RDGeom::Point3D pos3 = conf.getAtomPos(neighbor2->getIdx());

      // Calculate the vectors
      RDGeom::Point3D v1 = pos2 - pos1;
      RDGeom::Point3D v2 = pos3 - pos1;

      // Normalize the vectors
      v1.normalize();
      v2.normalize();

      // Calculate the angle in radians
      double angle = acos(v1.dotProduct(v2));

      // Convert to degrees
      double degrees = angle * 180.0 / M_PI;

      avgDegrees += degrees;
      n++;
    }
  }

  if (n >= 1) {
    avgDegrees /= n;
  }

  return avgDegrees;
}

void AssignHybridization(RDKit::RWMol &mol, RDKit::Conformer &conf) {
  for (auto atomIt = mol.beginAtoms(); atomIt != mol.endAtoms(); ++atomIt) {
    RDKit::Atom *atom = *atomIt;
    double avgDegrees = AverageBondAngle(mol, atom, conf);

    if (avgDegrees > 155.0) {
      atom->setHybridization(RDKit::Atom::HybridizationType::SP3);
    } else if (avgDegrees <= 155.0 && avgDegrees > 115.0) {
      atom->setHybridization(RDKit::Atom::HybridizationType::SP2);
    }
   
    // taken from AtomIsInRing
    auto ringInfo = atom->getOwningMol().getRingInfo();
    if (!ringInfo->isSssrOrBetter()) {
      RDKit::MolOps::findSSSR(atom->getOwningMol());
    }
    auto inRing = ringInfo->numAtomRings(atom->getIdx()) != 0;

    // special case for imines
    // RDKit::PeriodicTable::getTable()->getAtomicNumber("N")) {
    if (atom->getAtomicNum() == 7 && atom->getNumExplicitHs() == 1 &&
        atom->getDegree() == 2 && avgDegrees > 109.5) {
      atom->setHybridization(RDKit::Atom::HybridizationType::SP2);
    } else if (atom->getAtomicNum() == 7 && atom->getDegree() == 2 && inRing) { // azete
      atom->setHybridization(RDKit::Atom::HybridizationType::SP2);
    } // pass 1
    //


    // log out all conditions
    std::cout << "conds: " 
              << atom->getAtomicNum() << " "
              << atom->getNumExplicitHs() << " "
              << atom->getDegree() << " "
              << inRing << " "
              << hybridizationToString(atom->getHybridization()) << " "
              << std::endl;
  }
}

void WIP(RDKit::RWMol &mol, RDKit::Conformer &conf) {

  // RDKit::Atom *atom;
  for (auto atomIt = mol.beginAtoms(); atomIt != mol.endAtoms(); ++atomIt) {
    // atom = *atomIt;

    // Loop through bonds to find the initial maxbond
    for (auto bondIt = mol.getAtomBonds(*atomIt); bondIt.first != bondIt.second;
         ++bondIt.first) {
    }
  }
}

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
  //
  //

  std::cout << "Number of bonds w/ _struct_conns " << mol.getNumBonds()
            << std::endl;

  // print average bond angle for first 10 atoms
  // int i = 0;
  // for (const auto &atom : mol.atoms()) {
  //   if (i >= 10) {
  //     break;
  //   }
  //   std::cout << "Atom: " << atom->getIdx()
  //             << " Avg. Bond Angle: " << AverageBondAngle(mol, atom, *conf)
  //             << std::endl;
  //   i++;
  // }

  double angleTotal = 0.0;
  for (const auto &atom : mol.atoms()) {
    angleTotal += AverageBondAngle(mol, atom, *conf);
  }
  std::cout << "ANGLE TOTAL: " << angleTotal << std::endl;

  // if (!atom->getOwningMol().getRingInfo()->isSssrOrBetter()) {
  //   MolOps::findSSSR(atom->getOwningMol());
  // }

  RDKit::MolOps::setHybridization(mol);
  AssignHybridization(mol, *conf);

  return 0;
}
