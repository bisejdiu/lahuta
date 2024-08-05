#include "convert.hpp"
#include "GraphMol/MonomerInfo.h"
#include "GraphMol/RDKitBase.h"

#define ITER_GEMMI_ATOMS(st, atom)                                             \
  for (const Model &model : st.models)                                         \
    for (const Chain &chain : model.chains)                                    \
      for (const Residue &res : chain.residues)                                \
        for (const gemmi::Atom &atom : res.atoms)

using namespace gemmi;
using namespace RDKit;

// TODO: (@bis):
// 1. Residue names are truncated to 3 characters by gemmi. This may cause
// issues with some residues.
// 1. Cleanup the code and remove unnecessary comments
void gemmiStructureToRDKit(RWMol &mol, const Structure &st, Conformer &conf,
                           bool ign_h) {

  ITER_GEMMI_ATOMS(st, atom) {

    if (atom.element == Element("H") && ign_h) {
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

// TODO: there might be a better way to filter atoms from an RDKit molecule
RWMol rdMolFromRDKitMol(RWMol &mol, std::vector<int> &atomIndices) {
  Conformer conf = mol.getConformer();
  RWMol newMol;
  Conformer *newMolConf = new Conformer();

  for (auto atomIdx : atomIndices) {
    auto atom = mol.getAtomWithIdx(atomIdx);
    RDKit::Atom *newAtom = new RDKit::Atom(atom->getAtomicNum());
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

  return newMol;
}
