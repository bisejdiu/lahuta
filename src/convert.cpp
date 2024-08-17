#include <rdkit/GraphMol/MonomerInfo.h>
#include "convert.hpp"

#define ITER_GEMMI_ATOMS(st, atom)                                             \
  for (const Model &model : st.models)                                         \
    for (const Chain &chain : model.chains)                                    \
      for (const Residue &res : chain.residues)                                \
        for (const gemmi::Atom &atom : res.atoms)

using namespace gemmi;
using namespace RDKit;

void gemmiStructureToRDKit(RWMol &mol, const Structure &st, Conformer &conf,
                           bool ign_h) {

  ign_h = false; // FIX: ign_h=true is broken
  ITER_GEMMI_ATOMS(st, atom) {

    if (atom.element == Element("H") && ign_h) { // FIX: faster to swap the order
      continue;
    }

    int atomic_number = atom.element.atomic_number();
    RDKit::Atom *rAtom = new RDKit::Atom(atomic_number);
    rAtom->setFormalCharge(static_cast<int>(atom.charge));

    mol.addAtom(rAtom, true, true);

    auto pos = RDGeom::Point3D(atom.pos.x, atom.pos.y, atom.pos.z);
    conf.setAtomPos(atom.idx, pos);

    auto altLoc = (atom.altloc == '\0') ? "" : std::string(1, atom.altloc);
    AtomPDBResidueInfo atomInfo = {atom.name, atom.serial,         altLoc,
                                   res.name,  res.seqid.num.value, chain.name};

    atomInfo.setIsHeteroAtom(res.het_flag == 'H');
    atomInfo.setMonomerType(AtomMonomerInfo::PDBRESIDUE);

    auto *copy = static_cast<AtomMonomerInfo *>(atomInfo.copy());
    rAtom->setMonomerInfo(copy);
  }
}

RWMol rdMolFromRDKitMol(RWMol &mol, std::vector<int> &atomIndices) {
  Conformer conf = mol.getConformer();
  RWMol newMol;
  Conformer *newMolConf = new Conformer();

  for (auto atomIdx : atomIndices) {
    auto atom = mol.getAtomWithIdx(atomIdx);
    RDKit::Atom *newAtom = new RDKit::Atom(atom->getAtomicNum());
    // RDKit::Atom *newAtom = new RDKit::Atom(*atom); // FIX: check if this copies the monomer info
    newAtom->setFormalCharge(atom->getFormalCharge());
    newMol.addAtom(newAtom, true, true);
    auto pos = conf.getAtomPos(atom->getIdx());
    newMolConf->setAtomPos(newAtom->getIdx(), pos);

    auto *res = dynamic_cast<AtomPDBResidueInfo *>(atom->getMonomerInfo());
    AtomPDBResidueInfo atomInfo = {atom->getMonomerInfo()->getName(),
                                   static_cast<int>(atom->getIdx()),
                                   res->getAltLoc(),
                                   res->getResidueName(),
                                   res->getResidueNumber(),
                                   res->getChainId()};
    atomInfo.setIsHeteroAtom(res->getIsHeteroAtom());
    atomInfo.setMonomerType(AtomMonomerInfo::PDBRESIDUE);

    auto *copy = static_cast<AtomMonomerInfo *>(atomInfo.copy());
    newAtom->setMonomerInfo(copy);
  }

  newMol.addConformer(newMolConf, true);
  newMol.updatePropertyCache(false);

  return newMol;
}

RWMol rdMolFromRDKitMol(RWMol &mol, std::vector<int> &atomIndices, bool with_bonds) {
  if (!with_bonds) {
    return rdMolFromRDKitMol(mol, atomIndices);
  }
  Conformer conf = mol.getConformer();
  RWMol newMol;
  Conformer *newMolConf = new Conformer();

  for (auto atomIdx : atomIndices) {
    auto atom = mol.getAtomWithIdx(atomIdx);
    if (atom->getAtomicNum() == 1) {
      continue;
    }
    RDKit::Atom *newAtom = new RDKit::Atom(atom->getAtomicNum());
    newAtom->setFormalCharge(atom->getFormalCharge());
    newAtom->setIsAromatic(atom->getIsAromatic());
    newMol.addAtom(newAtom, true, true);
    auto pos = conf.getAtomPos(atom->getIdx());
    newMolConf->setAtomPos(newAtom->getIdx(), pos);

    auto *res = dynamic_cast<AtomPDBResidueInfo *>(atom->getMonomerInfo());
    AtomPDBResidueInfo atomInfo = {atom->getMonomerInfo()->getName(),
                                   static_cast<int>(atom->getIdx()),
                                   res->getAltLoc(),
                                   res->getResidueName(),
                                   res->getResidueNumber(),
                                   res->getChainId()};
    atomInfo.setIsHeteroAtom(res->getIsHeteroAtom());
    atomInfo.setMonomerType(AtomMonomerInfo::PDBRESIDUE);

    auto *copy = static_cast<AtomMonomerInfo *>(atomInfo.copy());
    newAtom->setMonomerInfo(copy);

    // std::cout << "Adding atom: " << newAtom->getIdx() << " " << newAtom->getSymbol() << std::endl;
    // // Add bonds to the newMol
    // for (auto bondIt = mol.getAtomBonds(atom); bondIt.first != bondIt.second;
    //      ++bondIt.first) {
    //   const RDKit::Bond *bond = mol[*bondIt.first];
    //   auto otherAtom = bond->getOtherAtom(atom);
    //   if (std::find(atomIndices.begin(), atomIndices.end(), otherAtom->getIdx()) != atomIndices.end()) {
    //     int aIx = std::distance(atomIndices.begin(), std::find(atomIndices.begin(), atomIndices.end(), atom->getIdx()));
    //     int bIx = std::distance(atomIndices.begin(), std::find(atomIndices.begin(), atomIndices.end(), otherAtom->getIdx()));
    //     if (newMol.getBondBetweenAtoms(aIx, bIx) == nullptr) {
    //       newMol.addBond(aIx, bIx, bond->getBondType());
    //     }
    //   }
    // }


  }

  for (auto atomIdx : atomIndices) {
    auto atom = mol.getAtomWithIdx(atomIdx);
    if (atom->getAtomicNum() == 1) {
      continue;
    }

    // std::cout << "Adding atom: " << atom->getIdx() << " " << atom->getSymbol() << std::endl;
    // Add bonds to the newMol
    for (auto bondIt = mol.getAtomBonds(atom); bondIt.first != bondIt.second;
         ++bondIt.first) {
      const RDKit::Bond *bond = mol[*bondIt.first];
      auto otherAtom = bond->getOtherAtom(atom);
      if (otherAtom->getAtomicNum() == 1) {
        continue;
      }
      if (std::find(atomIndices.begin(), atomIndices.end(), otherAtom->getIdx()) != atomIndices.end()) {
        int aIx = std::distance(atomIndices.begin(), std::find(atomIndices.begin(), atomIndices.end(), atom->getIdx()));
        int bIx = std::distance(atomIndices.begin(), std::find(atomIndices.begin(), atomIndices.end(), otherAtom->getIdx()));
        if (newMol.getBondBetweenAtoms(aIx, bIx) == nullptr) {
          newMol.addBond(aIx, bIx, bond->getBondType());
        }
      }
    }


  }

  newMol.addConformer(newMolConf, true);





  newMol.updatePropertyCache(false);
  return newMol;

}
