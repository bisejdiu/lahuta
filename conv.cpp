#include "conv.hpp"

using namespace gemmi;

RDKit::RWMol gemmiStructureToRDKit(Structure st, RDKit::Conformer &conf,
                                   bool ign_h) {

  RDKit::RWMol mol;

  int atom_index = 0;
  for (Model &model : st.models) {
    for (Chain &chain : model.chains) {
      // for (ResidueSpan& sub : chain.subchains()) {
      // if (const Entity* ent = st.get_entity_of(sub)) {
      for (Residue &res : chain.residues) {
        RDKit::AtomPDBResidueInfo res_info;
        res_info.setResidueName(res.name);
        res_info.setIsHeteroAtom(res.het_flag == 'H');
        for (const Atom &atom : res.atoms) {
          // std::cout << "Atom name: " << atom.name << std::endl;
          Element element = atom.element;
          if (element == Element("H") && ign_h) {
            continue;
          }
          // TODO: needs to be handled properly when ign_h is true
          if (element == Element("D")) {
            element = El::H;
          } else if (element == El::X) {
            element = El::X;
          }
          //
          int atomic_number = element.atomic_number();
          RDKit::Atom rdkit_atom(atomic_number);
          rdkit_atom.setFormalCharge((int)atom.charge);
          //
          auto *copy = (RDKit::AtomPDBResidueInfo *)res_info.copy();
          // copy->setName(res_name);
          // copy->setIsHeteroAtom((bool)res.het_flag);
          rdkit_atom.setMonomerInfo(copy);
          mol.addAtom(&rdkit_atom, true, false);

          conf.setAtomPos(atom_index,
                          RDGeom::Point3D(atom.pos.x, atom.pos.y, atom.pos.z));
          atom_index += 1;
        }
      }
    }
  }

  return mol;
}
