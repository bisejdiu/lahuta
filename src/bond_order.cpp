#include <GraphMol/RDKitBase.h>
#include "bond_order.hpp"
#include "ob/bond_utils.hpp"
#include "ob/kekulize.h"

using namespace RDKit;

using HybridizationType = RDKit::Atom::HybridizationType;
using SubStrMatches = std::vector<RDKit::MatchVectType>;

void PerceiveBondOrders(RDKit::RWMol &mol) {

  RDKit::Conformer &conf = mol.getConformer();

  // Set all atoms to sp hybridization (1 for openbabel)
  for (auto atomIt = mol.beginAtoms(); atomIt != mol.endAtoms(); ++atomIt) {
    auto atom = *atomIt;
    atom->setHybridization(HybridizationType::SP);
  }

  SubstructMatchParameters params;
  params.maxMatches = mol.getNumAtoms();
  // NOTE: SmartsMatching can take the number of threads as a parameter
  // params.numThreads = 1;
  RDKitSmartsMatch(mol, params);

  // Pass 1: Assign estimated hybridization based on average bond angles
  // taken from AtomIsInRing
  auto ringInfo = mol.getRingInfo();
  if (!ringInfo->isSssrOrBetter()) {
    RDKit::MolOps::findSSSR(mol);
  }
  for (auto atomIt = mol.beginAtoms(); atomIt != mol.endAtoms(); ++atomIt) {
    RDKit::Atom *atom = *atomIt;
    double avgDegrees = AverageBondAngle(atom);

    if (avgDegrees > 155.0) {
      atom->setHybridization(HybridizationType::SP);
    } else if (avgDegrees <= 155.0 && avgDegrees > 115.0) {
      atom->setHybridization(HybridizationType::SP2);
    }

    auto inRing = ringInfo->numAtomRings(atom->getIdx()) != 0;

    // special case for imines
    // RDKit::PeriodicTable::getTable()->getAtomicNumber("N")) {
    if (atom->getAtomicNum() == 7 && atom->getNumExplicitHs() == 1 &&
        atom->getDegree() == 2 && avgDegrees > 109.5) {
      atom->setHybridization(HybridizationType::SP2);
    } else if (atom->getAtomicNum() == 7 && atom->getDegree() == 2 &&
               inRing) { // azete
      atom->setHybridization(HybridizationType::SP2);
    } // pass 1
    //
  }

  // Pass 2: look for 5-memmbered rings with torsions <= 7.5 degrees
  //         and 6-membered rings with torsions <= 12 degrees
  //         (set all atoms with at least two bonds to sp2)

  if (!mol.getRingInfo()->isInitialized()) {
    RDKit::MolOps::findSSSR(mol);
  }

  double torsions = 0.0;
  auto rlist = mol.getRingInfo()->atomRings();
  for (const auto &ring : rlist) {
    if (ring.size() == 5) {
      std::vector<int> atom_indices;
      for (const auto &atomIdx : ring) {
        atom_indices.push_back(atomIdx);
      }

      RDKit::Atom *a0, *a1, *a2, *a3, *a4;
      a0 = mol.getAtomWithIdx(atom_indices[0]);
      a1 = mol.getAtomWithIdx(atom_indices[1]);
      a2 = mol.getAtomWithIdx(atom_indices[2]);
      a3 = mol.getAtomWithIdx(atom_indices[3]);
      a4 = mol.getAtomWithIdx(atom_indices[4]);

      // Get the coordinates of the atoms
      RDGeom::Point3D pos0 = conf.getAtomPos(a0->getIdx());
      RDGeom::Point3D pos1 = conf.getAtomPos(a1->getIdx());
      RDGeom::Point3D pos2 = conf.getAtomPos(a2->getIdx());
      RDGeom::Point3D pos3 = conf.getAtomPos(a3->getIdx());
      RDGeom::Point3D pos4 = conf.getAtomPos(a4->getIdx());

      torsions = (fabs(CalcTorsionAngle(pos0, pos1, pos2, pos3)) +
                  fabs(CalcTorsionAngle(pos1, pos2, pos3, pos4)) +
                  fabs(CalcTorsionAngle(pos2, pos3, pos4, pos0)) +
                  fabs(CalcTorsionAngle(pos3, pos4, pos0, pos1)) +
                  fabs(CalcTorsionAngle(pos4, pos0, pos1, pos2))) /
                 5.0;

      if (torsions <= 7.5) {
        for (const auto &atomIdx : ring) {
          RDKit::Atom *atom = mol.getAtomWithIdx(atomIdx);
          // if (atom->getDegree() == 2) {
          if (ob_explicit_valence(mol, atom) == 2) {
            atom->setHybridization(HybridizationType::SP2);
          }
        }
      }
    } else if (ring.size() == 6) {
      std::vector<int> atom_indices;
      for (const auto &atomIdx : ring) {
        atom_indices.push_back(atomIdx);
      }

      RDKit::Atom *a0, *a1, *a2, *a3, *a4, *a5;
      a0 = mol.getAtomWithIdx(atom_indices[0]);
      a1 = mol.getAtomWithIdx(atom_indices[1]);
      a2 = mol.getAtomWithIdx(atom_indices[2]);
      a3 = mol.getAtomWithIdx(atom_indices[3]);
      a4 = mol.getAtomWithIdx(atom_indices[4]);
      a5 = mol.getAtomWithIdx(atom_indices[5]);

      // Get the coordinates of the atoms
      RDGeom::Point3D pos0 = conf.getAtomPos(a0->getIdx());
      RDGeom::Point3D pos1 = conf.getAtomPos(a1->getIdx());
      RDGeom::Point3D pos2 = conf.getAtomPos(a2->getIdx());
      RDGeom::Point3D pos3 = conf.getAtomPos(a3->getIdx());
      RDGeom::Point3D pos4 = conf.getAtomPos(a4->getIdx());
      RDGeom::Point3D pos5 = conf.getAtomPos(a5->getIdx());

      torsions = (fabs(CalcTorsionAngle(pos0, pos1, pos2, pos3)) +
                  fabs(CalcTorsionAngle(pos1, pos2, pos3, pos4)) +
                  fabs(CalcTorsionAngle(pos2, pos3, pos4, pos5)) +
                  fabs(CalcTorsionAngle(pos3, pos4, pos5, pos0)) +
                  fabs(CalcTorsionAngle(pos4, pos5, pos0, pos1)) +
                  fabs(CalcTorsionAngle(pos5, pos0, pos1, pos2))) /
                 6.0;

      if (torsions <= 12.0) {
        for (const auto &atomIdx : ring) {
          RDKit::Atom *atom = mol.getAtomWithIdx(atomIdx);
          // if (atom->getDegree() == 2 || atom->getDegree() == 3) {
          if (ob_explicit_valence(mol, atom) == 2 ||
              ob_explicit_valence(mol, atom) == 3) {
            atom->setHybridization(HybridizationType::SP2);
          }
        }
      }
    }
  }

  // Pass 3: "Antialiasing" If an atom marked as sp hybrid isn't
  //          bonded to another or an sp2 hybrid isn't bonded
  //          to another (or terminal atoms in both cases)
  //          mark them to a lower hybridization for now

  // NOTE: Hybridization assignment for all atoms is done before function is
  // called
  for (auto atomIt = mol.beginAtoms(); atomIt != mol.endAtoms(); ++atomIt) {
    RDKit::Atom *atom = *atomIt;
    if (atom->getHybridization() == HybridizationType::SP ||
        atom->getHybridization() == HybridizationType::SP2) {
      bool openNbr = false;
      ROMol::ADJ_ITER nbrIdx, endNbr;
      boost::tie(nbrIdx, endNbr) = mol.getAtomNeighbors(atom);
      for (; nbrIdx != endNbr; ++nbrIdx) {
        const RDKit::Atom *neighbor = mol.getAtomWithIdx(*nbrIdx);
        if (neighbor->getHybridization() < HybridizationType::SP3 ||
            ob_explicit_valence(mol, neighbor) == 1) {
          openNbr = true;
          break;
        }
      }
      if (!openNbr && atom->getHybridization() == HybridizationType::SP2) {
        atom->setHybridization(HybridizationType::SP3);
      } else if (!openNbr &&
                 atom->getHybridization() == HybridizationType::SP) {
        atom->setHybridization(HybridizationType::SP2);
      }
    }
  }
  // pass 3

  // Pass 4: Check for known functional group patterns and assign bonds
  //         to the canonical form
  //      Currently we have explicit code to do this, but a "bond typer"
  //      is in progress to make it simpler to test and debug.

  OBBondTypeAssignment(mol);

  std::string carbo("[#8D1;!-][#6](*)(*)");
  RDKit::ROMol *pattern = RDKit::SmartsToMol(carbo);

  std::vector<RDKit::MatchVectType> matchList;
  RDKit::SubstructMatch(mol, *pattern, matchList);

  for (const auto &match : matchList) {
    RDKit::Atom *a1 = mol.getAtomWithIdx(match[0].second);
    RDKit::Atom *a2 = mol.getAtomWithIdx(match[1].second);

    double avgDegrees = AverageBondAngle(a1);
    auto a1Pos = conf.getAtomPos(a1->getIdx());
    auto a2Pos = conf.getAtomPos(a2->getIdx());
    RDGeom::Point3D v1 = a2Pos - a1Pos;
    double dist1 = v1.length();

    // carbonyl geometries ?
    if (avgDegrees > 115.0 && avgDegrees < 150.0 && dist1 < 1.28) {
      bool hasDoubleBond = false;
      for (auto bondIt = mol.getAtomBonds(a1); bondIt.first != bondIt.second;
           ++bondIt.first) {
        const RDKit::Bond *bond = mol[*bondIt.first];
        if (bond->getBondType() == Bond::BondType::DOUBLE) {
          hasDoubleBond = true;
          break;
        }
      }
      if (!hasDoubleBond) {
        RDKit::Bond *bond = mol.getBondBetweenAtoms(a1->getIdx(), a2->getIdx());
        if (bond) {
          bond->setBondType(Bond::BondType::DOUBLE);
        }
      }
    }

  } // Carbonyl oxygen

  // thione C=S
  std::string thione("[#16D1][#6](*)(*)");
  RDKit::ROMol *pattern2 = RDKit::SmartsToMol(thione);

  std::vector<RDKit::MatchVectType> matchList2;
  RDKit::SubstructMatch(mol, *pattern2, matchList2);

  for (const auto &match : matchList2) {
    RDKit::Atom *a1 = mol.getAtomWithIdx(match[0].second);
    RDKit::Atom *a2 = mol.getAtomWithIdx(match[1].second);

    double avgDegrees = AverageBondAngle(a1);
    auto a1Pos = conf.getAtomPos(a1->getIdx());
    auto a2Pos = conf.getAtomPos(a2->getIdx());
    RDGeom::Point3D v1 = a2Pos - a1Pos;
    double dist1 = v1.length();

    // thione geometries ?
    if (avgDegrees > 115.0 && avgDegrees < 150.0 && dist1 < 1.72) {
      bool hasDoubleBond = false;
      for (auto bondIt = mol.getAtomBonds(a1); bondIt.first != bondIt.second;
           ++bondIt.first) {
        const RDKit::Bond *bond = mol[*bondIt.first];
        if (bond->getBondType() == Bond::BondType::DOUBLE) {
          hasDoubleBond = true;
          break;
        }
      }
      if (!hasDoubleBond) {
        RDKit::Bond *bond = mol.getBondBetweenAtoms(a1->getIdx(), a2->getIdx());
        if (bond) {
          bond->setBondType(Bond::BondType::DOUBLE);
        }
      }
    }
  } // thione

  // Isocyanate N=C=O or Isothiocyanate
  std::string isocyanate("[#8,#16;D1][#6D2][#7D2]");
  RDKit::ROMol *pattern3 = RDKit::SmartsToMol(isocyanate);

  std::vector<RDKit::MatchVectType> matchList3;
  RDKit::SubstructMatch(mol, *pattern3, matchList3);

  for (const auto &match : matchList3) {
    RDKit::Atom *a1 = mol.getAtomWithIdx(match[0].second);
    RDKit::Atom *a2 = mol.getAtomWithIdx(match[1].second);
    RDKit::Atom *a3 = mol.getAtomWithIdx(match[2].second);

    double avgDegrees = AverageBondAngle(a2);
    auto a1Pos = conf.getAtomPos(a1->getIdx());
    auto a2Pos = conf.getAtomPos(a2->getIdx());
    auto a3Pos = conf.getAtomPos(a3->getIdx());
    RDGeom::Point3D v1 = a2Pos - a1Pos;
    RDGeom::Point3D v2 = a3Pos - a2Pos;
    double dist1 = v1.length();
    double dist2 = v2.length();

    // isocyanate geometry or Isothiocyanate geometry ?
    double dist1OK;
    if (a1->getAtomicNum() == 8)
      auto dist1OK = dist1 < 1.28;
    else
      auto dist1OK = dist1 < 1.72;

    if (avgDegrees > 150.0 && dist1OK && dist2 < 1.34) {
      RDKit::Bond *bond1 = mol.getBondBetweenAtoms(a1->getIdx(), a2->getIdx());
      RDKit::Bond *bond2 = mol.getBondBetweenAtoms(a2->getIdx(), a3->getIdx());

      if (!bond1 || !bond2)
        continue;
      bond1->setBondType(Bond::BondType::DOUBLE);
      bond2->setBondType(Bond::BondType::DOUBLE);
    }
  } // Isocyanate

  // oxime C=S
  std::string oxime("[#6D3][#7D2][#8D2]");
  RDKit::ROMol *pattern4 = RDKit::SmartsToMol(oxime);

  std::vector<RDKit::MatchVectType> matchList4;
  RDKit::SubstructMatch(mol, *pattern4, matchList4);

  for (const auto &match : matchList4) {
    RDKit::Atom *a1 = mol.getAtomWithIdx(match[0].second);
    RDKit::Atom *a2 = mol.getAtomWithIdx(match[1].second);

    double avgDegrees = AverageBondAngle(a1);
    auto a1Pos = conf.getAtomPos(a1->getIdx());
    auto a2Pos = conf.getAtomPos(a2->getIdx());
    RDGeom::Point3D v1 = a2Pos - a1Pos;
    double dist1 = v1.length();

    // thione geometries ?
    if (avgDegrees > 110.0 && avgDegrees < 150.0 && dist1 < 1.4) {
      bool hasDoubleBond = false;
      for (auto bondIt = mol.getAtomBonds(a1); bondIt.first != bondIt.second;
           ++bondIt.first) {
        const RDKit::Bond *bond = mol[*bondIt.first];
        if (bond->getBondType() == Bond::BondType::DOUBLE) {
          hasDoubleBond = true;
          break;
        }
      }
      if (!hasDoubleBond) {
        RDKit::Bond *bond = mol.getBondBetweenAtoms(a1->getIdx(), a2->getIdx());
        if (bond) {
          bond->setBondType(Bond::BondType::DOUBLE);
        }
      }
    }
  } // oxime

  // oxido-n+ (e.g., pyridine-N-oxide)
  std::string oxidopyr("[#8D1][#7D3r6]");
  RDKit::ROMol *pattern5 = RDKit::SmartsToMol(oxidopyr);

  std::vector<RDKit::MatchVectType> matchList5;
  RDKit::SubstructMatch(mol, *pattern5, matchList5);

  for (const auto &match : matchList5) {
    RDKit::Atom *a1 = mol.getAtomWithIdx(match[0].second);
    RDKit::Atom *a2 = mol.getAtomWithIdx(match[1].second);

    double avgDegrees = AverageBondAngle(a1);
    auto a1Pos = conf.getAtomPos(a1->getIdx());
    auto a2Pos = conf.getAtomPos(a2->getIdx());
    RDGeom::Point3D v1 = a2Pos - a1Pos;
    double dist1 = v1.length();

    if (avgDegrees > 110.0 && avgDegrees < 150.0 && dist1 < 1.35) {
      a1->setFormalCharge(-1); // oxygen
      a2->setFormalCharge(+1); // nitrogen
    }
  } // oxido-n+

  // Pass 5: Check for aromatic rings and assign bonds as appropriate
  // This is just a quick and dirty approximation that marks everything
  //  as potentially aromatic

  // This doesn't work perfectly, but it's pretty decent.
  //  Need to have a list of SMARTS patterns for common rings
  //  which would "break ties" on complicated multi-ring systems
  // (Most of the current problems lie in the interface with the
  //   Kekulize code anyway, not in marking everything as potentially aromatic)
  bool needs_kekulization = false;
  bool typed = false;
  unsigned int loopSize;
  for (const auto &ring : rlist) {
    typed = false;
    loopSize = ring.size();
    if (loopSize == 5 || loopSize == 6 || loopSize == 7) {
      for (const auto &atomIdx : ring) {
        RDKit::Atom *atom = mol.getAtomWithIdx(atomIdx);
        bool hasBondOfOrder2 = false;
        bool hasBondOfOrder3 = false;
        for (auto bondIt = mol.getAtomBonds(atom);
             bondIt.first != bondIt.second; ++bondIt.first) {
          const RDKit::Bond *bond = mol[*bondIt.first];
          if (bond->getBondType() == Bond::BondType::DOUBLE) {
            hasBondOfOrder2 = true;
          } else if (bond->getBondType() == Bond::BondType::TRIPLE) {
            hasBondOfOrder3 = true;
          }
        }
        if (hasBondOfOrder2 || hasBondOfOrder3 ||
            atom->getHybridization() != HybridizationType::SP2) {
          typed = true;
          break;
        }
      }

      if (!typed) {
        for (unsigned int loop = 0; loop < loopSize; ++loop) {
          RDKit::Bond *bond =
              mol.getBondBetweenAtoms(ring[loop], ring[(loop + 1) % loopSize]);
          if (bond) {
            bond->setIsAromatic(true);
            needs_kekulization = true;
          }
        }
      }
    }
  }

  if (needs_kekulization) {
    for (auto bondIt = mol.beginBonds(); bondIt != mol.endBonds(); ++bondIt) {
      RDKit::Bond *bond = *bondIt;
      if (bond->getIsAromatic()) {
        bond->getBeginAtom()->setIsAromatic(true);
        bond->getEndAtom()->setIsAromatic(true);
      }
    }
    bool ok = OBKekulize(&mol);
  }
}
