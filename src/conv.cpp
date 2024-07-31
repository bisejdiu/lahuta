#include "conv.hpp"
#include "GraphMol/MonomerInfo.h"
#include "GraphMol/RDKitBase.h"

#define ITER_GEMMI_ATOMS(st, atom)                                             \
  for (const Model &model : st.models)                                         \
    for (const Chain &chain : model.chains)                                    \
      for (const Residue &res : chain.residues)                                \
        for (const Atom &atom : res.atoms)

using namespace gemmi;

// TODO: (@bis):
// 1. Residue names are truncated to 3 characters by gemmi. This may cause
// issues with some residues.
// 1. Cleanup the code and remove unnecessary comments
void gemmiStructureToRDKit(RDKit::RWMol &mol, const Structure &st,
                           RDKit::Conformer &conf, bool ign_h) {

  // RDKit::RWMol mol;

  int aIx = 0;
  ITER_GEMMI_ATOMS(st, atom) {
    // FIX: should Hydrogens be added or ignored?
    // if (atom.altloc != '\0') {
    //   continue;
    // }
    // std::cout << "ALTLOC: " << atom.altloc << std::endl;
    Element element = atom.element;
    if (element == Element("H") && ign_h) {
      continue;
    }
    // FIX: Check how RDKit handles this (e.g. PDBAtomFromSymbol)
    // TODO: needs to be handled properly when ign_h is true
    // if (element == Element("D")) {
    //   element = El::H;
    // } else if (element == El::X) {
    //   element = El::X;
    // }
    //

    int atomic_number = element.atomic_number();
    RDKit::Atom *rAtom = new RDKit::Atom(atomic_number);
    rAtom->setFormalCharge(static_cast<int>(atom.charge));

    mol.addAtom(rAtom, true, true);

    auto pos = RDGeom::Point3D(atom.pos.x, atom.pos.y, atom.pos.z);
    conf.setAtomPos(aIx, pos);

    const auto &atomName = atom.name;
    // std::string altLoc(1, atom.altloc) means that altLoc is a string of
    // length 1, with the first character being atom.altloc
    // This is done to convert char to string
    std::string altLoc =
        (atom.altloc == '\0') ? "" : std::string(1, atom.altloc);
    // std::cout << "res.name" << res.name << std::endl;
    RDKit::AtomPDBResidueInfo atomInfo = {
        atom.name, atom.serial,         altLoc,
        res.name,  res.seqid.num.value, chain.name};
    bool het_flag = res.het_flag == 'H';
    atomInfo.setIsHeteroAtom(het_flag);
    atomInfo.setMonomerType(RDKit::AtomMonomerInfo::PDBRESIDUE);

    RDKit::AtomMonomerInfo *copy =
        static_cast<RDKit::AtomMonomerInfo *>(atomInfo.copy());
    rAtom->setMonomerInfo(copy);

    aIx += 1;
  }
}

// FIX: Not efficient, because we have to re-find the atoms in the new molecule
RDKit::RWMol rdMolFromRDKitMol(RDKit::RWMol &mol,
                               std::vector<int> &atomIndices) {
  RDKit::Conformer conf = mol.getConformer();
  RDKit::RWMol newMol;
  RDKit::Conformer *newMolConf = new RDKit::Conformer();
  // for (auto atomIt = mol.beginAtoms(); atomIt != mol.endAtoms(); ++atomIt) {
  //   auto atom = *atomIt;
  //   auto it = std::find(atomIndices.begin(), atomIndices.end(),
  //   atom->getIdx()); if (it != atomIndices.end()) {
  //     RDKit::Atom *newAtom = new RDKit::Atom(atom->getAtomicNum());
  //     newAtom->setFormalCharge(atom->getFormalCharge());
  //     newMol.addAtom(newAtom, true, true);
  //     auto pos = conf.getAtomPos(atom->getIdx());
  //     newMolConf->setAtomPos(newAtom->getIdx(), pos);
  //
  //     RDKit::AtomPDBResidueInfo *res = dynamic_cast<RDKit::AtomPDBResidueInfo
  //     *>(atom->getMonomerInfo()); std::string atomName =
  //     atom->getMonomerInfo()->getName(); std::string altLoc =
  //     res->getAltLoc(); std::string resName = res->getResidueName(); int
  //     resSeq = res->getResidueNumber(); std::string chainId =
  //     res->getChainId(); RDKit::AtomPDBResidueInfo atomInfo = {
  //         atomName, static_cast<int>(atom->getIdx()), altLoc, resName,
  //         resSeq, chainId};
  //     atomInfo.setIsHeteroAtom(res->getIsHeteroAtom());
  //     atomInfo.setMonomerType(RDKit::AtomMonomerInfo::PDBRESIDUE);
  //
  //     RDKit::AtomMonomerInfo *copy =
  //         static_cast<RDKit::AtomMonomerInfo *>(atomInfo.copy());
  //     newAtom->setMonomerInfo(copy);
  //
  //   }
  // }

  for (auto atomIdx : atomIndices) {
    auto atom = mol.getAtomWithIdx(atomIdx);
    RDKit::Atom *newAtom = new RDKit::Atom(atom->getAtomicNum());
    newAtom->setFormalCharge(atom->getFormalCharge());
    newMol.addAtom(newAtom, true, true);
    auto pos = conf.getAtomPos(atom->getIdx());
    newMolConf->setAtomPos(newAtom->getIdx(), pos);

    RDKit::AtomPDBResidueInfo *res =
        dynamic_cast<RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    std::string atomName = atom->getMonomerInfo()->getName();
    std::string altLoc = res->getAltLoc();
    std::string resName = res->getResidueName();
    int resSeq = res->getResidueNumber();
    std::string chainId = res->getChainId();
    RDKit::AtomPDBResidueInfo atomInfo = {
        atomName, static_cast<int>(atom->getIdx()), altLoc, resName, resSeq,
        chainId};
    atomInfo.setIsHeteroAtom(res->getIsHeteroAtom());
    atomInfo.setMonomerType(RDKit::AtomMonomerInfo::PDBRESIDUE);

    RDKit::AtomMonomerInfo *copy =
        static_cast<RDKit::AtomMonomerInfo *>(atomInfo.copy());
    newAtom->setMonomerInfo(copy);
  }

  newMol.addConformer(newMolConf, true);

  return newMol;
}
