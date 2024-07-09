#include "Geometry/point.h"
#include "GraphMol/Atom.h"
#include "GraphMol/Bond.h"
#include "GraphMol/DetermineBonds/DetermineBonds.h"
#include "GraphMol/MonomerInfo.h"
#include "GraphMol/PeriodicTable.h"
#include "GraphMol/RDKitBase.h"
#include "GraphMol/RWMol.h"

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

#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "elements.h"
#include "kekulize.h"

#include <iostream>

using namespace gemmi;
using namespace RDKit;

using HybridizationType = RDKit::Atom::HybridizationType;
using SubStrMatches = std::vector<RDKit::MatchVectType>;

int GetExpVal(const RDKit::Atom *atom) {
  double valence = 0;
  for (auto bondIt = atom->getOwningMol().getAtomBonds(atom);
       bondIt.first != bondIt.second; ++bondIt.first) {
    const RDKit::Bond *bond = atom->getOwningMol()[*bondIt.first];
    valence += bond->getBondTypeAsDouble();
  }
  return valence;
}

constexpr std::pair<const char *, HybridizationType> smartsList[] = {
    {"[D4]", HybridizationType::SP3},
    {"[D5]", HybridizationType::SP3D},
    {"[D6]", HybridizationType::SP3D2},
    {"[C]", HybridizationType::SP3},
    {"[c,$(C=*)]", HybridizationType::SP2},
    {"[$(C#*),$(C(=*)=*)]", HybridizationType::SP},

    {"[N]", HybridizationType::SP3},
    {"[n,$(N=*),$(N[#6,#7,#8]=,:,#*)]", HybridizationType::SP2},
    {"[ND1,ND2,ND3]a", HybridizationType::SP2},
    {"[$(N#*),$([ND2](=*)=*)]", HybridizationType::SP},

    {"[O]", HybridizationType::SP3},
    {"[o,$(O=*),$(O[#6,#7,#8]=,:*)]", HybridizationType::SP2},
    {"[$([#8D1][#6][#8D1])]", HybridizationType::SP2},
    {"[$(O#*)]", HybridizationType::SP},

    {"[S]", RDKit::Atom::HybridizationType::SP3},
    {"[#16;s,$([SD1]=*)]", RDKit::Atom::HybridizationType::SP2},
    {"[SD6]", RDKit::Atom::HybridizationType::SP3D2},
};

SubStrMatches performSubstructMatch(RDKit::ROMol &mol, RDKit::ROMol &pattern,
                                    SubstructMatchParameters &params) {
  SubStrMatches matchList;
  matchList = RDKit::SubstructMatch(mol, pattern, params);

  if (matchList.size() == params.maxMatches) {
    SubstructMatchParameters largeParams = params;
    largeParams.maxMatches = mol.getNumAtoms();
    matchList = RDKit::SubstructMatch(mol, pattern, largeParams);
  }

  return matchList;
}

void RDKitSmartsMatch(RDKit::ROMol &mol, SubstructMatchParameters &params) {
  // Precompute patterns statically
  static std::array<RDKit::ROMol *, std::size(smartsList)> patterns = [] {
    std::array<RDKit::ROMol *, std::size(smartsList)> temp{};
    for (size_t i = 0; i < std::size(smartsList); ++i) {
      temp[i] = RDKit::SmartsToMol(smartsList[i].first);
    }
    return temp;
  }();

  for (size_t i = 0; i < std::size(smartsList); ++i) {
    const auto &[smarts, hybridType] = smartsList[i];
    RDKit::ROMol *pattern = patterns[i];

    // FIX: supplying the param list leads to performance issues
    // SubStrMatches matchList = performSubstructMatch(mol, *pattern, params);
    SubStrMatches matchList;
    RDKit::SubstructMatch(mol, *pattern, matchList);

    for (const auto &match : matchList) {
      auto atom = mol.getAtomWithIdx(match[0].second);
      atom->setHybridization(hybridType);
    }
  }
}

using SmartsPatternPair = std::pair<const char *, std::vector<int>>;
const std::vector<SmartsPatternPair> bondSmarts = {
    {"[x2,x3]1[#6]([#7D3]2)[#6][#6][#6]2[x2,x3][#6]([#7D3]3)[#6][#6][#6]3["
     "x2,"
     "x3][#6]([#7D3]4)[#6][#6][#6]4[x2,x3][#6]([#7D3]5)[#6][#6][#6]51",
     {0,  1,  2,  1,  2,  1,  1,  3,  1,  3,  4,  2,  4,  5,  1,  5,  2,
      1,  5,  6,  2,  6,  7,  1,  7,  8,  2,  7,  9,  1,  9,  10, 2,  10,
      11, 1,  11, 8,  1,  11, 12, 2,  12, 13, 1,  13, 14, 1,  13, 15, 2,
      15, 16, 1,  16, 17, 2,  17, 14, 1,  17, 18, 1,  18, 19, 2,  19, 20,
      1,  19, 21, 1,  21, 22, 2,  22, 23, 1,  23, 20, 2}},
    {"[x2,x3]1[#6]([#7D3]2)[#6][#6][#6]2[x2,x3][#6]([#7]3)[#6][#6][#6]3[x2,"
     "x3][#6]([#7D3]4)[#6][#6][#6]4[x2,x3][#6]([#7]5)[#6][#6][#6]51",
     {0,  1,  2,  1,  2,  1,  1,  3,  1,  3,  4,  2,  4,  5,  1,  5,  2,
      1,  5,  6,  2,  6,  7,  1,  7,  8,  2,  7,  9,  1,  9,  10, 2,  10,
      11, 1,  11, 8,  1,  11, 12, 2,  12, 13, 1,  13, 14, 1,  13, 15, 2,
      15, 16, 1,  16, 17, 2,  17, 14, 1,  17, 18, 1,  18, 19, 2,  19, 20,
      1,  19, 21, 1,  21, 22, 2,  22, 23, 1,  23, 20, 2}},
    {"[x2,x3]1[#6]([#7]2)[#6][#6][#6]2[x2,x3][#6]([#7]3)[#6][#6][#6]3[x2,"
     "x3]["
     "#6]([#7]4)[#6][#6][#6]4[x2,x3][#6]([#7]5)[#6][#6][#6]51",
     {0,  1,  2,  1,  2,  1,  1,  3,  1,  3,  4,  2,  4,  5,  1,  5,  2,
      1,  5,  6,  2,  6,  7,  1,  7,  8,  2,  7,  9,  1,  9,  10, 2,  10,
      11, 1,  11, 8,  1,  11, 12, 2,  12, 13, 1,  13, 14, 1,  13, 15, 2,
      15, 16, 1,  16, 17, 2,  17, 14, 1,  17, 18, 1,  18, 19, 2,  19, 20,
      1,  19, 21, 1,  21, 22, 2,  22, 23, 1,  23, 20, 2}},
    {"[#7D2][#7D2^1][#7D1]", {0, 1, 2, 1, 2, 2}},
    {"[#8D1][#7D3^2]([#8D1])*", {0, 1, 2, 1, 2, 2, 1, 3, 1}},
    {"[#16D4]([#8D1])([#8D1])([*!#8])([*!#8])",
     {0, 1, 2, 0, 2, 2, 0, 3, 1, 0, 4, 1}},
    {"[#16D4]([#8D1])([#8D1])([#8-,#8D1])([#8-,#8D1])",
     {0, 1, 2, 0, 2, 2, 0, 3, 1, 0, 4, 1}},
    {"[#16D4]([#16D1])([#8D1])([#8-,#8])([#8-,#8])",
     {0, 1, 2, 0, 2, 2, 0, 3, 1, 0, 4, 1}},
    {"[#16D3]([#8D1])([*!#8])([*!#8])", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#16D3]([#8D1])([#8D1-])([#8D1-])", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#16D3]([#8D1])([#8])([#8])", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#16D2]([#8D1])([#16D1])", {0, 1, 2, 0, 2, 2}},
    {"[#16D2]([#8D1])([*!#8])", {0, 1, 2, 0, 2, 1}},
    {"[#16D2]([#8D1])([#8D1])", {0, 1, 2, 0, 2, 2}},
    {"[#15D3]([#8D1])([#8D1])([#8D2])", {0, 1, 2, 0, 2, 2, 0, 3, 1}},
    {"[#7D2]([#8D1])([#1])", {0, 1, 2, 0, 2, 1}},
    {"[#15D4]([#8D1])(*)(*)(*)", {0, 1, 2, 0, 2, 1, 0, 3, 1, 0, 4, 1}},
    {"[#6D3^2]([#8D1])([#8])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#8D1][#6D2^1][#8D1]", {0, 1, 2, 1, 2, 2}},
    {"[#6D3^2]([#8D1;!-])([#7])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#34D3^2]([#8D1])([#8])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#6D3^2]([#8D1])([#16])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#6D3^2]([#16D1])([#16])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[CD3^2]([#16D1])([N])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    // Do not seem to work correctly with RDKit:
    // {"[#6^2][#6D2^1][#6^2]", {0, 1, 2, 1, 2, 2}},
    // {"[#6^2][#6D2^1][#8D1]", {0, 1, 2, 1, 2, 2}},
    {"[#6D3^2;!R]([#7D1H0;!R])([#7;!R])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#6D3^2;!R]([#7D2H1;!R])([#7;!R])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#6D3^2;!R]([#7D3H2;!R])([#7;!R])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#6D3^2;!R]([#1,#6])([#1,#6])[#7D3^2;!R]([#1])[#6]",
     {0, 1, 1, 0, 2, 1, 0, 3, 2, 3, 4, 1, 3, 5, 1}},
};

void OBBondTypeAssignment(RDKit::ROMol &mol) {

  // Precompute patterns statically
  static std::vector<RDKit::ROMol *> patterns = [] {
    std::vector<RDKit::ROMol *> temp;
    temp.resize(bondSmarts.size());
    for (size_t i = 0; i < bondSmarts.size(); ++i) {
      temp[i] = RDKit::SmartsToMol(bondSmarts[i].first);
    }
    return temp;
  }();

  for (size_t i = 0; i < bondSmarts.size(); ++i) {
    const auto &[smarts, bondVector] = bondSmarts[i];
    RDKit::ROMol *pattern = patterns[i];

    SubstructMatchParameters params;
    params.maxMatches = 100000000;
    std::vector<RDKit::MatchVectType> matchList;
    matchList = RDKit::SubstructMatch(mol, *pattern, params);

    // debug:
    // if (matchList.size() != 0) {
    //   std::cout << "Match Success with: " << smarts << " " <<
    //   matchList.size()
    //             << std::endl;
    // }

    for (const auto &match : matchList) {
      for (auto j = 0; j < bondVector.size(); j += 3) {
        auto bond = mol.getBondBetweenAtoms(match[bondVector[j]].second,
                                            match[bondVector[j + 1]].second);
        if (bond) {
          // FIX: there is an implicit conversion from int to BondType
          bond->setBondType(Bond::BondType(bondVector[j + 2]));
        }
      }
    }
  }
}

double AverageBondAngle(const RDKit::Atom *atom) {
  double avgDegrees = 0.0;
  int n = 0;

  const RDKit::ROMol &mol = atom->getOwningMol();
  const RDKit::Conformer &conf = mol.getConformer();

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

// FIX: RDKit likely has a function for this
//
/*!  This function calculates the torsion angle of three vectors, represented
  by four points A--B--C--D, i.e. B and C are vertexes, but none of A--B,
  B--C, and C--D are colinear.  A "torsion angle" is the amount of "twist"
  or torsion needed around the B--C axis to bring A--B into the same plane
  as B--C--D.  The torsion is measured by "looking down" the vector B--C so
  that B is superimposed on C, then noting how far you'd have to rotate
  A--B to superimpose A over D.  Angles are + in the anti-clockwise
  direction.  The operation is symmetrical in that if you reverse the image
  (look from C to B and rotate D over A), you get the same answer.
*/

double CalcTorsionAngle(const RDGeom::Point3D &a, const RDGeom::Point3D &b,
                        const RDGeom::Point3D &c, const RDGeom::Point3D &d) {

  double torsion;
  RDGeom::Point3D b1, b2, b3, c1, c2, c3;

  b1 = a - b;
  b2 = b - c;
  b3 = c - d;

  double rb2 = sqrt(b2.dotProduct(b2));

  RDGeom::Point3D b2xb3 = b2.crossProduct(b3);
  RDGeom::Point3D b1xb2 = b1.crossProduct(b2);
  torsion = -atan2(rb2 * b1.dotProduct(b2xb3), b1xb2.dotProduct(b2xb3));

  return (torsion * 180.0 / M_PI);
}

void PerceiveBondOrders(RDKit::RWMol &mol) {

  RDKit::Conformer &conf = mol.getConformer();

  // Pass 1: Assign estimated hybridization based on average bond angles
  for (auto atomIt = mol.beginAtoms(); atomIt != mol.endAtoms(); ++atomIt) {
    RDKit::Atom *atom = *atomIt;
    double avgDegrees = AverageBondAngle(atom);

    if (avgDegrees > 155.0) {
      atom->setHybridization(HybridizationType::SP);
    } else if (avgDegrees <= 155.0 && avgDegrees > 115.0) {
      atom->setHybridization(HybridizationType::SP2);
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
    std::cout << "SSSR not initialized" << std::endl;
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
          if (GetExplicitValenceFromAtomBonds(mol, atom) == 2) {
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
          if (GetExplicitValenceFromAtomBonds(mol, atom) == 2 ||
              GetExplicitValenceFromAtomBonds(mol, atom) == 3) {
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
            GetExplicitValenceFromAtomBonds(mol, neighbor) == 1) {
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
  //

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

  // NOTE: the valences do not differ, but RDKit atoms & bonds need updating,
  // before getting the correct valence for (const auto &atom: mol.atoms()) {
  //   auto valence = GetExplicitValenceFromAtomBonds(mol, atom);
  //   if (atom->getExplicitValence() != valence) {
  //     std::cout << "V Diff: Atom: " << atom->getIdx() << " Valence: " <<
  //     atom->getExplicitValence() << " " << valence << std::endl;
  //   }
  // }

  // Clean up incorrect bonds
  // start = std::chrono::high_resolution_clock::now();
  // CleanUpMolecule(mol, *conf);
  // end = std::chrono::high_resolution_clock::now();
  // elapsed = end - start;
  // std::cout << "Time to clean up molecule: " << elapsed.count() * 1000 << "
  // ms"
  //           << std::endl;

  // NOTE: should we clean up bonds before or after adding the struct conns?
  mol.updatePropertyCache(false);
  // Clean up incorrect bonds
  start = std::chrono::high_resolution_clock::now();
  CleanUpMolecule(mol, *conf);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to clean up molecule: " << elapsed.count() * 1000 << "ms "
            << std::endl;

  for (gemmi::Connection &conn : st.connections) {
    // FIX: Need to iterate over all models
    gemmi::Atom *a1 = st.first_model().find_cra(conn.partner1).atom;
    gemmi::Atom *a2 = st.first_model().find_cra(conn.partner2).atom;

    if (mol.getBondBetweenAtoms(a1->serial - 1, a2->serial - 1) == nullptr) {
      mol.addBond((unsigned int)a1->serial - 1, (unsigned int)a2->serial - 1,
                  RDKit::Bond::BondType::SINGLE);
      continue;
    } 
  }

  std::cout << "Number of bonds w/ _struct_conns " << mol.getNumBonds()
            << std::endl;

  // get ROMol from RWMol using cast
  // auto _mol = static_cast<RDKit::ROMol*>(new RDKit::RWMol(mol));

  for (auto atomIt = mol.beginAtoms(); atomIt != mol.endAtoms(); ++atomIt) {
    auto atom = *atomIt;
    atom->setHybridization(HybridizationType::SP);
  }
  start = std::chrono::high_resolution_clock::now();
  SubstructMatchParameters params;
  // params.maxMatches = 500; // 00000;
  // params.numThreads = 1;
  RDKitSmartsMatch(mol, params);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to match SMARTS pattern: " << elapsed.count() * 1000
            << " ms" << std::endl;

  start = std::chrono::high_resolution_clock::now();
  PerceiveBondOrders(mol);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to perceive bond orders: " << elapsed.count() * 1000
            << " ms" << std::endl;

  // for (auto atom : mol.atoms()) {
  //   auto resi = atom->getMonomerInfo();
  //   AtomPDBResidueInfo *residue = dynamic_cast<AtomPDBResidueInfo *>(resi);
  //   std::cout << "Atom Info: " << atom->getIdx() << " " << atom->getSymbol()
  //             << " " << resi->getMonomerType() << " " << resi->getName() << "
  //             "
  //             << residue->getResidueName() << " " <<
  //             residue->getIsHeteroAtom()
  //             << " " << std::endl;
  // }

  std::string bondOrders = "";
  for (auto bondIt = mol.beginBonds(); bondIt != mol.endBonds(); ++bondIt) {
    RDKit::Bond *bond = *bondIt;
    bondOrders += std::to_string(bond->getBeginAtomIdx()) + " " +
                  std::to_string(bond->getEndAtomIdx()) + " " +
                  std::to_string(bond->getBondType()) + "\n";
  }
  std::cout << "FINAL RESULT: " << std::endl;
  std::cout << bondOrders << std::endl;

  return 0;
}
