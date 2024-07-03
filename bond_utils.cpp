#include "Geometry/point.h"
#include "GraphMol/Atom.h"
#include "GraphMol/Bond.h"
#include "GraphMol/RDKitBase.h"
#include "GraphMol/RWMol.h"

#include "elements.h"
#include "bond_utils.hpp"

double computeLengthSq(const RDKit::Conformer &conf, const RDKit::Bond *bond) {
  auto a1Pos = conf.getAtomPos(bond->getBeginAtomIdx());
  auto a2Pos = conf.getAtomPos(bond->getEndAtomIdx());
  // compute the length of the bond
  double d2x = (a1Pos.x - a2Pos.x);
  double d2y = (a1Pos.y - a2Pos.y);
  double d2z = (a1Pos.z - a2Pos.z);

  return d2x * d2x + d2y * d2y + d2z * d2z;
}

unsigned int GetExplicitValenceFromAtomBonds(const RDKit::RWMol &mol,
                                             const RDKit::Atom *atom) {
  unsigned int valence = 0;
  for (auto bond : boost::make_iterator_range(mol.getAtomBonds(atom))) {
    valence++; // when this is called, all bonds are single bonds
  }
  return valence;
}

bool isAngleLessThan45Degrees(const RDKit::RWMol &mol,
                              const RDKit::Conformer &conf,
                              const RDKit::Atom *atom1,
                              const RDKit::Atom *atom2,
                              const RDKit::Atom *atom3) {

  RDGeom::Point3D v1, v2;
  RDGeom::Point3D p1 = conf.getAtomPos(atom1->getIdx());
  RDGeom::Point3D p2 = conf.getAtomPos(atom2->getIdx());
  RDGeom::Point3D p3 = conf.getAtomPos(atom3->getIdx());

  // return isAngleLessThan45(p1, p2, p3);
  v1 = p1 - p2;
  v2 = p3 - p2;

  double lenSq1 = v1.lengthSq();
  double lenSq2 = v2.lengthSq();

  if (lenSq1 < 1e-8 || lenSq2 < 1e-8) {
    return false; // One of the vectors is too small
  }

  double dotProd = v1.dotProduct(v2);

  // Check if the angle is less than 45 degrees using the dot product directly
  // cos(45 degrees) = sqrt(2) / 2, so we need to check if
  // dotProd > sqrt(lenSq1 * lenSq2) * sqrt(2) / 2
  // This can be simplified to:
  // 2 * dotProd^2 > lenSq1 * lenSq2
  return dotProd > 0 && 2.0 * dotProd * dotProd > lenSq1 * lenSq2;
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
  // double angle = vectorAngle(v1, v2);

  return angle;
}

double ComputeSmallestBondAngle(const RDKit::RWMol &mol, const RDKit::Conformer &conf,
                       const RDKit::Atom *atom) {
  double minDegrees = 360.0;

  for (auto bondIt = mol.getAtomBonds(atom); bondIt.first != bondIt.second;
       ++bondIt.first) {
    auto bond1 = mol[*bondIt.first];
    auto atom1 = bond1->getOtherAtom(atom);
    for (auto bondIt2 = mol.getAtomBonds(atom); bondIt2.first != bondIt2.second;
         ++bondIt2.first) {
      if (bondIt.first == bondIt2.first)
        continue;
      auto bond2 = mol[*bondIt2.first];
      auto atom2 = bond2->getOtherAtom(atom);

      double angle = computeAngle(mol, conf, atom1, atom, atom2);
      if (angle < minDegrees) {
        minDegrees = angle;
      }
    }
  }

  return minDegrees;
}
bool SmallestBondAngle(const RDKit::RWMol &mol, const RDKit::Conformer &conf,
                       const RDKit::Atom *atom) {
  double minDegrees = 360.0;

  for (auto bondIt = mol.getAtomBonds(atom); bondIt.first != bondIt.second;
       ++bondIt.first) {
    auto bond1 = mol[*bondIt.first];
    auto atom1 = bond1->getOtherAtom(atom);
    for (auto bondIt2 = mol.getAtomBonds(atom); bondIt2.first != bondIt2.second;
         ++bondIt2.first) {
      if (bondIt.first == bondIt2.first)
        continue;
      auto bond2 = mol[*bondIt2.first];
      auto atom2 = bond2->getOtherAtom(atom);

      if (isAngleLessThan45Degrees(mol, conf, atom1, atom, atom2)) {
        return true;
      }
    }
  }

  return false;
}

// NOTE: The most time consuming part is removing bonds from the molecule, which
// is out of our control. RDKit seems to provide a batch bond removal utility,
// which may provide a performance boost, be we'd nned to refactor the function
// to use it (e.g. keeping track of the explicit valence and remembering what
// bonds have been removed)
void CleanUpMolecule(RDKit::RWMol &mol, RDKit::Conformer &conf) {
  RDKit::Atom *atom;
  RDKit::Bond *maxbond, *bond;
  double maxlength;
  int valCount;
  bool changed;

  for (auto atomIt = mol.beginAtoms(); atomIt != mol.endAtoms(); ++atomIt) {
    atom = *atomIt;
    while (GetExplicitValenceFromAtomBonds(mol, atom) >
               OBElements::GetMaxBonds(atom->getAtomicNum()) ||
           // SmallestBondAngle(mol, conf, atom)) {
           ComputeSmallestBondAngle(mol, conf, atom) < 45.0) {

      maxbond = nullptr;

      // Loop through bonds to find the initial maxbond
      for (auto bondIt = mol.getAtomBonds(atom); bondIt.first != bondIt.second;
           ++bondIt.first) {
        bond = mol[*bondIt.first];
        maxbond = bond;
        break;
      }

      if (maxbond == nullptr) {
        break; // No bonds to process
      }

      // // Delete bonds between hydrogens when over max valence
      if (atom->getAtomicNum() == 1) { // Hydrogen
        // std::cout << "Deleting Hydrogen bond: " << atom->getIdx() <<
        // std::endl;
        changed = false;
        for (auto bondIt = mol.getAtomBonds(atom);
             bondIt.first != bondIt.second; ++bondIt.first) {
          bond = mol[*bondIt.first];
          if (bond->getOtherAtom(atom)->getAtomicNum() ==
              1) { // Neighboring Hydrogen
            mol.removeBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
            changed = true;
            break;
          }
        }
        if (changed) {
          continue; // Reevaluate
        }
      }

      auto maxlength = computeLengthSq(conf, maxbond);
      int i = 0;
      for (auto bondIt = mol.getAtomBonds(atom); bondIt.first != bondIt.second;
           ++bondIt.first) {
        if (i == 0) {
          i++;
          continue;
        }
        bond = mol[*bondIt.first];
        auto length = computeLengthSq(conf, bond);
        if (length > maxlength) {
          maxbond = bond;
          maxlength = length;
        }
      }

      // auto atom1_name = atom->getMonomerInfo()->getName();
      // auto atom2_name =
      // maxbond->getOtherAtom(atom)->getMonomerInfo()->getName(); std::cout <<
      // "Deleting bond: " << atom->getIdx() << " - " <<
      // maxbond->getOtherAtom(atom)->getIdx() << " " << atom1_name << " - " <<
      // atom2_name << std::endl;

      mol.removeBond(maxbond->getBeginAtomIdx(), maxbond->getEndAtomIdx());
    }
  }
}
