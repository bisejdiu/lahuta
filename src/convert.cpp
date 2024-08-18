#include <rdkit/GraphMol/MonomerInfo.h>
#include <unordered_map>
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

// FIX: It seems RDKit automatically copies bonds or coordinates when adding 
// atoms directly. Need to investigate this. 
RWMol filter_atoms(RWMol &mol, std::vector<int> &atomIndices) {
  RWMol newMol;

  for (auto atomIdx : atomIndices) {
    auto atom = mol.getAtomWithIdx(atomIdx);
    RDKit::Atom *newAtom = new RDKit::Atom(*atom); 
    newMol.addAtom(newAtom, true, true);
  }

  newMol.updatePropertyCache(false);

  return newMol;
}

RWMol filter_atom_conf(RWMol &mol, std::vector<int> &atomIndices) {
  Conformer conf = mol.getConformer();
  RWMol newMol;
  Conformer *newMolConf = new Conformer();

  for (auto atomIdx : atomIndices) {
    auto atom = mol.getAtomWithIdx(atomIdx);
    RDKit::Atom *newAtom = new RDKit::Atom(*atom); 
    newMol.addAtom(newAtom, true, true);
    auto pos = conf.getAtomPos(atom->getIdx());
    newMolConf->setAtomPos(newAtom->getIdx(), pos);
  }

  newMol.addConformer(newMolConf, true);
  newMol.updatePropertyCache(false);

  return newMol;
}

RWMol filter_with_atom_data(RWMol &mol, std::vector<int> &atomIndices) {
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

RWMol rdMolFromRDKitMol(RWMol &mol, std::vector<int> &indices) {

  RWMol new_mol;
  Conformer conf = mol.getConformer();
  Conformer *new_conf = new Conformer();

  std::unordered_set<unsigned int> _ixs; 
  auto get_or_add_atom = [&](const RDKit::Atom *at, unsigned int &_ix) {

    if (_ixs.find(_ix) == _ixs.end()) {
      RDKit::Atom *atom = new RDKit::Atom(*at);
      new_mol.addAtom(atom, true, true);
      _ixs.insert(_ix);

      auto pos = conf.getAtomPos(atom->getIdx());
      new_conf->setAtomPos(atom->getIdx(), pos);
    }
    return _ix++; 
  };

  unsigned int _ix = 0;
  for (auto idx: indices) {
    auto at = mol.getAtomWithIdx(idx);
    auto a_ix = get_or_add_atom(at, _ix);

    for (auto bondIt = mol.getAtomBonds(at); bondIt.first != bondIt.second; ++bondIt.first) {
      const RDKit::Bond *bond = mol[*bondIt.first];
      auto oat = bond->getOtherAtom(at);
      if (std::find(indices.begin(), indices.end(), oat->getIdx()) != indices.end()) {
        auto b_ix = get_or_add_atom(oat, _ix);
        if (new_mol.getBondBetweenAtoms(a_ix, b_ix) == nullptr) {
          new_mol.addBond(a_ix, b_ix, bond->getBondType());
        }
      }
    }
  }

  new_mol.addConformer(new_conf, true);
  new_mol.updatePropertyCache(false);
  return new_mol;

}

