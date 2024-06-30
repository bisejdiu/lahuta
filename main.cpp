#include "GraphMol/RDKitBase.h"
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
#include <gemmi/model.hpp> // for Model
// #include <gemmi/assembly.hpp>
#include <gemmi/mmread_gz.hpp> // for read_structure_gz
// #include <gemmi/neighbor.hpp>

#include "bonds.hpp"
#include "conv.hpp"
#include "nsgrid.hpp"

#include <chrono>
#include <gemmi/metadata.hpp>

#include "elements.h"

using namespace gemmi;

double vectorAngle(const RDGeom::Point3D &v1, const RDGeom::Point3D &v2) {
  double dp;

  dp = v1.dotProduct(v2) / (v1.length() * v2.length());
  if (dp < -0.999999)
    dp = -0.999999;

  if (dp > 0.999999)
    dp = 0.999999;

  return acos(dp) * 180.0 / M_PI;
}

double computeAngle(const RDKit::RWMol &mol, const RDKit::Conformer &conf, const RDKit::Atom *atom1, const RDKit::Atom *atom2,
                    RDKit::Atom *atom3) {

  RDGeom::Point3D v1, v2;
  RDGeom::Point3D p1 = conf.getAtomPos(atom1->getIdx());
  RDGeom::Point3D p2 = conf.getAtomPos(atom2->getIdx());
  RDGeom::Point3D p3 = conf.getAtomPos(atom3->getIdx());

  v1 = p1 - p2;
  v2 = p3 - p2;

  if (v1.length() < 0.0001 || v2.length() < 0.0001) {
    return 0.0;
  } 

  double angle = v1.angleTo(v2) * 180.0 / M_PI;

  return angle;
}


double SmallestBondAngle(const RDKit::RWMol &mol, const RDKit::Conformer &conf, const RDKit::Atom *atom) {
  double minDegrees = 360.0;

  // RDKit::ROMol::OBOND_ITER_PAIR atomBonds = atom->getOwningMol().getAtomBonds(atom);
  // while (atomBonds.first != atomBonds.second) {
  //   unsigned int bondIdx = atom->getOwningMol().getTopology()[*atomBonds.first]->getIdx();
  //   RDKit::Bond *bond = mol.getBondWithIdx(bondIdx);
  //   RDKit::Atom *neighbor = bond->getOtherAtom(atom);
  //
  //   for (auto bond2 : boost::make_iterator_range(mol.getAtomBonds(atom))) {
  //     if (bond->getIdx() == bond2) continue;
  //     auto atom2 = mol[bond2]->getOtherAtom(atom);
  //
  //     double angle = computeAngle(mol, conf, neighbor, atom, atom2);
  //     if (angle < minDegrees) {
  //       minDegrees = angle;
  //     }
  //   }
  //
  // }
  // return 0;


  for (auto bond1 : boost::make_iterator_range(mol.getAtomBonds(atom))) {
    auto atom1 = mol[bond1]->getOtherAtom(atom);
    for (auto bond2 : boost::make_iterator_range(mol.getAtomBonds(atom))) {
      if (bond1 == bond2) continue;
      auto atom2 = mol[bond2]->getOtherAtom(atom);

      double angle = computeAngle(mol, conf, atom1, atom, atom2);
      if (angle < minDegrees) {
        minDegrees = angle;
      }
    }
  }

  return minDegrees;
}


// double SmallestBondAngle(const RDKit::RWMol &mol, RDKit::Conformer &conf, const RDKit::Atom *atom) {
//   double minDegrees = 360.0;
//
//   for (auto bond1 : boost::make_iterator_range(mol.getAtomBonds(atom))) {
//     auto atom1 = mol.getBondWithIdx(bond1)->getOtherAtom(atom);
//     for (auto bond2 : boost::make_iterator_range(mol.getAtomBonds(atom))) {
//       if (bond1 == bond2) continue;
//       auto atom2 = mol.getBondWithIdx(bond2)->getOtherAtom(atom);
//
//       double angle = computeAngle(mol, conf, atom1, atom, atom2);
//       if (angle < minDegrees) {
//         minDegrees = angle;
//       }
//     }
//   }
//
//   return minDegrees;
// }


void CleanUpMolecule(RDKit::RWMol &mol, RDKit::Conformer &conf) {
  RDKit::Atom *atom;
  RDKit::Bond *maxbond, *bond;
  double maxlength;
  int valCount;
  bool changed;

  for (auto atomIt = mol.beginAtoms(); atomIt != mol.endAtoms(); ++atomIt) {
    atom = *atomIt;
    // while (atom->getTotalDegree() > RDKit::PeriodicTable::getTable()->getNouterElecs(atom->getAtomicNum())
    //        || SmallestBondAngle(mol, conf, atom) < 45.0) {
    while (atom->getTotalDegree() > OBElements::GetMaxBonds(atom->getAtomicNum())) {
      maxbond = nullptr;
      valCount = 0;

      // Loop through bonds to find the initial maxbond
      for (auto bondIt = mol.getAtomBonds(atom); bondIt.first != bondIt.second; ++bondIt.first) {
        bond = mol[*bondIt.first];
        maxbond = bond;
        valCount++;
      }

      if (valCount == 0) {
        break; // No bonds to process
      }

      // Delete bonds between hydrogens when over max valence
      if (atom->getAtomicNum() == 1) { // Hydrogen
        changed = false;
        for (auto bondIt = mol.getAtomBonds(atom); bondIt.first != bondIt.second; ++bondIt.first) {
          bond = mol[*bondIt.first];
          if (bond->getOtherAtom(atom)->getAtomicNum() == 1) { // Neighboring Hydrogen
            mol.removeBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
            changed = true;
            break;
          }
        }
        if (changed) {
          continue; // Reevaluate
        }
      }

      // Find the bond with the maximum length
      auto a1Pos = conf.getAtomPos(maxbond->getBeginAtomIdx());
      auto a2Pos = conf.getAtomPos(maxbond->getEndAtomIdx());
      // compute the length of the bond
        // d2 = SQUARE(begin->GetX() - end->GetX());
        // d2 += SQUARE(begin->GetY() - end->GetY());
        // d2 += SQUARE(begin->GetZ() - end->GetZ());
        // return(sqrt(d2));
      auto d2 = (a1Pos.x - a2Pos.x) * (a1Pos.x - a2Pos.x);
      d2 += (a1Pos.y - a2Pos.y) * (a1Pos.y - a2Pos.y);
      d2 += (a1Pos.z - a2Pos.z) * (a1Pos.z - a2Pos.z);
      double length = sqrt(d2);

      
      for (auto bondIt = mol.getAtomBonds(atom); bondIt.first != bondIt.second; ++bondIt.first) {
        bond = mol[*bondIt.first];
        a1Pos = conf.getAtomPos(bond->getBeginAtomIdx());
        a2Pos = conf.getAtomPos(bond->getEndAtomIdx());
        // compute the length of the bond
        // d2 = SQUARE(begin->GetX() - end->GetX());
        // d2 += SQUARE(begin->GetY() - end->GetY());
        // d2 += SQUARE(begin->GetZ() - end->GetZ());
        // return(sqrt(d2));
        d2 = (a1Pos.x - a2Pos.x) * (a1Pos.x - a2Pos.x);
        d2 += (a1Pos.y - a2Pos.y) * (a1Pos.y - a2Pos.y);
        d2 += (a1Pos.z - a2Pos.z) * (a1Pos.z - a2Pos.z);
        double length = sqrt(d2);

        if (length > maxlength) {
          maxbond = bond;
          maxlength = length;
        }
      }

      // Delete the bond with the longest length
      mol.removeBond(maxbond->getBeginAtomIdx(), maxbond->getEndAtomIdx());
    }
  }
}


// double SmallestBondAngle(RDKit::RWMol &mol, RDKit::Conformer &conf, RDKit::Atom *atom) {
//   double degrees, minDegrees;
//   auto bonds = mol.getAtomBonds(atom); // returns type: OBOND_ITER_PAIR
//   minDegrees = 360.0;
//
//   boost::detail::out_edge_iter first = bonds.first;
//   boost::detail::out_edge_iter second = bonds.second;
//
//   for (; first != second; ++first) {
//     RDKit::Bond *bond = mol[*first];
//
//     RDKit::Atom *neighbor = bond->getOtherAtom(atom);
//
//   }
//
//   return minDegrees;
// }

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


  // std::map<std::string, std::vector<std::string>> atoms;
  // gemmi::cif::Block& block = gemmi::cif::read(file_name).sole_block();
  // auto atoms = block.find_mmcif_category("_atom_site.");
  // time rdKit conversion

  start = std::chrono::high_resolution_clock::now();
  RDKit::Conformer *conf = new RDKit::Conformer();
  RDKit::RWMol mol = gemmiStructureToRDKit(st, *conf, false);
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
  Model &model = st.first_model();
  // time to find bonds
  start = std::chrono::high_resolution_clock::now();
  findBondsDeconstructed(st, model, mol, 5.0, 1.3);
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

  double cutoff = 5.001;

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
  perceiveBonds(mol, results, 0.45);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to find bonds using NSGrid: " << elapsed.count() * 1000 << " ms"
            << std::endl;
  std::cout << "Number of bonds found using NSGrid: " << mol.getNumBonds() << std::endl;

  // for (gemmi::Connection& conn: st.connections) {
  //   // std::cout << "_struc_conn: " << conn.partner1.str() << " " << conn.partner2.str() << std::endl;
  //   // gemmi::Model::find_cra(conn.partner1);
  //   // FIX: Need to iterate over all models
  //   Atom* a1 = st.first_model().find_cra(conn.partner1).atom;
  //   Atom* a2 = st.first_model().find_cra(conn.partner2).atom;
  //
  //   if (mol.getBondBetweenAtoms(a1->serial, a2->serial) != nullptr) {
  //     continue;
  //   }
  //   mol.addBond(a1->serial, a2->serial, RDKit::Bond::BondType::SINGLE);
  // }

  // for (const auto &atom: mol.atoms()) {
  //   // std::cout << "Atom name: " << atom->getSymbol() << " "
  //   //           << atom->getIdx() << std::endl;
  //   double val = SmallestBondAngle(mol, *conf, atom);
  //   while (atom->getExplicitValence() > OBElements::GetMaxBonds(atom->getAtomicNum()) || val < 90.0) {
  //
  //   }
  // }
  //

  CleanUpMolecule(mol, *conf);

  std::cout << "Number of bonds w/ _struct_conns " << mol.getNumBonds() << std::endl;
  mol.updatePropertyCache(false);
  std::cout << "Number of bonds w/ _struct_conns " << mol.getNumBonds() << std::endl;
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
