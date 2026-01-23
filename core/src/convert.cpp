#include <rdkit/GraphMol/MonomerInfo.h>

#include "convert.hpp"
#include "logging/logging.hpp"

namespace lahuta {

void add_atom_to_mol(RDKit::RWMol &mol, RDKit::Conformer &conf, const gemmi::Atom &atom, const gemmi::Residue &res, const gemmi::Chain &chain) {
  unsigned int atomic_number  = atom.element.atomic_number();
  RDKit::Atom *rdkit_atom = new RDKit::Atom(atomic_number);
  rdkit_atom->setFormalCharge(static_cast<int>(atom.charge));

  mol.addAtom(rdkit_atom, true, true);

  auto pos = RDGeom::Point3D(atom.pos.x, atom.pos.y, atom.pos.z);
  conf.setAtomPos(atom.idx, pos);

  auto alt_loc = (atom.altloc == '\0') ? "" : std::string(1, atom.altloc);
  // FIX: num is OptionalNum (not guaranteed to have a value)
  auto *info = new RDKit::AtomPDBResidueInfo(atom.name, atom.serial, alt_loc, res.name, res.seqid.num.value, chain.name);


  info->setResidueIndex(res.idx);
  info->setIsHeteroAtom(res.het_flag == 'H');
  info->setMonomerType(RDKit::AtomMonomerInfo::PDBRESIDUE);

  rdkit_atom->setMonomerInfo(info);
}

void create_RDKit_repr(RDKit::RWMol &mol, const gemmi::Structure &st, RDKit::Conformer &conf) {
  if (st.models.empty()) {
    Logger::get_logger()->warn("No models found in the structure");
    return;
  }

  if (st.models.size() > 1) {
    Logger::get_logger()->warn("Structure contains {} models, only processing the first one ({})",
                              st.models.size(), st.models[0].name);
  }

  const gemmi::Model &model = st.models[0];

  for (const gemmi::Chain &chain : model.chains) {
    for (const gemmi::Residue &res : chain.residues) {
      for (const gemmi::Atom &atom : res.atoms) {
        add_atom_to_mol(mol, conf, atom, res, chain);
      }
    }
  }
}

std::shared_ptr<RDKit::RWMol> create_RDKit(const gemmi::Structure &st) {
  auto mol = std::make_shared<RDKit::RWMol>();
  RDKit::Conformer *conformer = new RDKit::Conformer();
  create_RDKit_repr(*mol, st, *conformer);
  mol->updatePropertyCache(false);
  mol->addConformer(conformer, true);
  return mol;
}

void IR_to_RWMol(RDKit::RWMol &mol, const IR &ir) {
  RDKit::Conformer *conf = new RDKit::Conformer();
  for (size_t i = 0; i < ir.atom_indices.size(); ++i) {
    RDKit::Atom *atom = new RDKit::Atom(ir.atomic_numbers[i]);
    atom->setMonomerInfo(
        new RDKit::AtomPDBResidueInfo(ir.atom_names[i], -1, "", ir.resnames[i], ir.resids[i], ir.chainlabels[i]));
    mol.addAtom(atom, true, true);
    conf->setAtomPos(i, {ir.positions[i][0], ir.positions[i][1], ir.positions[i][2]});
  }
  mol.addConformer(conf, true);
  mol.updatePropertyCache(false);
}

} // namespace lahuta
