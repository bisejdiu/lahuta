#include "convert.hpp"
#include <rdkit/GraphMol/MonomerInfo.h>
#include <unordered_map>

// clang-format off
#define ITER_GEMMI_ATOMS(st, model, chain, res, atom)               \
  for (const gemmi::Model &model : st.models)                       \
    for (const gemmi::Chain &chain : model.chains)                  \
      for (const gemmi::Residue &res : chain.residues)              \
        for (const gemmi::Atom &atom : res.atoms)
// clang-format on

using namespace gemmi;
using namespace RDKit;

namespace lahuta {

// FIX: We should build the Topology from here
// FIX: Converters should in fact be loaders
void create_RDKit_repr(RWMol &mol, const Structure &st, Conformer &conf, bool ign_h) {

  ign_h = false; // FIX: ign_h=true is broken
  bool is_first_model = true;
  std::string curr_model = "";
  ITER_GEMMI_ATOMS(st, model, chain, res, atom) {

    // temporarily only processing the first model
    if (is_first_model) {
      is_first_model = false;
      curr_model = model.name; 
    } else if (curr_model != model.name) {
      break;
    }

    if (ign_h && atom.element == Element("H")) { // FIX: use `is_hydrogen` (elem.hpp) instead
      continue;
    }

    int atomic_number = atom.element.atomic_number();
    RDKit::Atom *rAtom = new RDKit::Atom(atomic_number);
    rAtom->setFormalCharge(static_cast<int>(atom.charge));

    mol.addAtom(rAtom, true, true);

    auto pos = RDGeom::Point3D(atom.pos.x, atom.pos.y, atom.pos.z);
    conf.setAtomPos(atom.idx, pos);

    auto altLoc = (atom.altloc == '\0') ? "" : std::string(1, atom.altloc);
    // FIX: num is OptionalNum (not guaranteed to have a value)
    AtomPDBResidueInfo atomInfo = {atom.name, atom.serial, altLoc, res.name, res.seqid.num.value, chain.name};

    atomInfo.setResidueIndex(res.idx);
    atomInfo.setIsHeteroAtom(res.het_flag == 'H');
    atomInfo.setMonomerType(AtomMonomerInfo::PDBRESIDUE);

    auto *copy = static_cast<AtomMonomerInfo *>(atomInfo.copy());
    rAtom->setMonomerInfo(copy);
  }
}

std::shared_ptr<RDKit::RWMol> create_RDKit(const Structure &st) {
  auto mol = std::make_shared<RDKit::RWMol>();
  RDKit::Conformer *conformer = new RDKit::Conformer();
  create_RDKit_repr(*mol, st, *conformer, false);
  mol->updatePropertyCache(false);
  mol->addConformer(conformer, true);
  return mol;
}

void IR_to_RWMol(RWMol &mol, const IR &ir) {
  Conformer *conf = new Conformer();
  for (size_t i = 0; i < ir.atom_indices.size(); ++i) {
    RDKit::Atom *atom = new RDKit::Atom(ir.atomic_numbers[i]);
    atom->setMonomerInfo(
        new AtomPDBResidueInfo(ir.atom_names[i], -1, "", ir.resnames[i], ir.resids[i], ir.chainlabels[i]));
    mol.addAtom(atom, true, true);
    conf->setAtomPos(i, {ir.positions[i][0], ir.positions[i][1], ir.positions[i][2]});
  }
  mol.addConformer(conf, true);
  mol.updatePropertyCache(false);
}

RWMol filter_atoms(RWMol &mol, std::vector<int> &indices) {
  RWMol new_mol;

  for (auto atomIdx : indices) {
    auto atom = mol.getAtomWithIdx(atomIdx);
    RDKit::Atom *new_atom = new RDKit::Atom(*atom);
    new_mol.addAtom(new_atom, true, true);
  }
  new_mol.updatePropertyCache(false);

  return new_mol;
}

RWMol filter_with_conf(RWMol &mol, std::vector<int> &indices) {
  Conformer conf = mol.getConformer();
  RWMol new_mol;
  Conformer *new_conf = new Conformer();

  for (auto atomIdx : indices) {
    auto atom = mol.getAtomWithIdx(atomIdx);
    RDKit::Atom *newAtom = new RDKit::Atom(*atom);
    new_mol.addAtom(newAtom, true, true);
    auto pos = conf.getAtomPos(atom->getIdx());
    new_conf->setAtomPos(newAtom->getIdx(), pos);
  }
  new_mol.addConformer(new_conf, true);
  new_mol.updatePropertyCache(false);

  return new_mol;
}

RWMol filter_with_bonds(const RWMol &mol, const std::vector<int> &indices) {
  RWMol new_mol;
  const Conformer &conf = mol.getConformer();
  auto *new_conf = new Conformer();

  // for faster index lookup
  std::unordered_map<int, int> old_to_new_index;

  // sanity check
  auto max_idx = std::max_element(indices.begin(), indices.end());
  if (max_idx == indices.end() || *max_idx >= mol.getNumAtoms()) {
    new_mol.addConformer(new_conf, true);
    new_mol.updatePropertyCache(false);
    return new_mol;
  }

  for (size_t i = 0; i < indices.size(); ++i) {
    int idx = indices[i];
    const RDKit::Atom *atom = mol.getAtomWithIdx(idx);
    int new_idx = new_mol.addAtom(new RDKit::Atom(*atom), false, true);
    old_to_new_index[idx] = new_idx;
    new_conf->setAtomPos(new_idx, conf.getAtomPos(idx));
  }

  for (int old_idx : indices) {
    for (const auto &bond : mol.atomBonds(mol.getAtomWithIdx(old_idx))) {
      int begin_idx = bond->getBeginAtomIdx();
      int end_idx = bond->getEndAtomIdx();

      auto begin_it = old_to_new_index.find(begin_idx);
      auto end_it = old_to_new_index.find(end_idx);

      if (begin_it != old_to_new_index.end() && end_it != old_to_new_index.end()) {
        int new_begin_idx = begin_it->second;
        int new_end_idx = end_it->second;
        if (new_begin_idx > new_end_idx) std::swap(new_begin_idx, new_end_idx);

        if (!new_mol.getBondBetweenAtoms(new_begin_idx, new_end_idx)) {
          new_mol.addBond(new_begin_idx, new_end_idx, bond->getBondType());
        }
      }
    }
  }

  new_mol.addConformer(new_conf, true);
  new_mol.updatePropertyCache(false);
  return new_mol;
}

} // namespace lahuta
