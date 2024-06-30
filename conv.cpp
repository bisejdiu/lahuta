// #include <gemmi/mmcif.hpp>
// #include <gemmi/atox.hpp>    // for string_to_int
// #include <gemmi/mmcif_impl.hpp> // for set_cell_from_mmcif
// // #include <gemmi/enumstr.hpp> // for entity_type_from_string, polymer_type_from_string
// #include <gemmi/polyheur.hpp>  // for restore_full_ccd_codes

#include "conv.hpp"

using namespace gemmi;

// void copy_int(const cif::Table::Row& row, int n, int& dest) {
//   if (row.has2(n))
//     dest = cif::as_int(row[n]);
// }
// void copy_double(const cif::Table::Row& row, int n, double& dest) {
//   if (row.has2(n))
//     dest = cif::as_number(row[n]);
// }
// void copy_string(const cif::Table::Row& row, int n, std::string& dest) {
//   if (row.has2(n))
//     dest = cif::as_string(row[n]);
// }
//
// enum Type { Covale=0, Disulf, Hydrog, MetalC, Unknown };
//
// std::string connection_type_to_string(Type t) {
//   switch (t) {
//     case Covale: return "covale";
//     case Disulf: return "disulf";
//     case Hydrog: return "hydrog";
//     case MetalC: return "metalc";
//     default: return "unknown";
//   }
// }
//
// Type connection_type_from_string(std::string& t) {
//   for (int i = 0; i != Connection::Unknown; ++i)
//     if (connection_type_to_string(Type(i)) == t)
//       return Type(i);
//   return Unknown;
// }
//
// void read_connectivity(cif::Block& block, Structure& st) {
//   enum {
//     kId=0, kConnTypeId=1,
//     kAuthAsymId=2/*-3*/,  kLabelAsymId=4/*-5*/, kLabelCompId=6/*-7*/,
//     kLabelAtomId=8/*-9*/, kLabelAltId=10/*-11*/,
//     kAuthSeqId=12/*-13*/, kLabelSeqId=14/*-15*/, kInsCode=16/*-17*/,
//     kSym1=18, kSym2=19, kDistValue=20, kLinkId=21
//   };
//   // label_ identifiers are not sufficient for HOH:
//   // waters have null label_seq_id so we need auth_seq_id+icode.
//   // And since we need auth_seq_id, we also use auth_asym_id for consistency.
//   // Unless only label_*_id are available.
//   for (const auto row : block.find("_struct_conn.", {
//         "id", "conn_type_id",                                   // 0-1
//         "?ptnr1_auth_asym_id", "?ptnr2_auth_asym_id",           // 2-3
//         "?ptnr1_label_asym_id", "?ptnr2_label_asym_id",         // 4-5
//         "ptnr1_label_comp_id", "ptnr2_label_comp_id",           // 6-7
//         "ptnr1_label_atom_id", "ptnr2_label_atom_id",           // 8-9
//         "?pdbx_ptnr1_label_alt_id", "?pdbx_ptnr2_label_alt_id", // 10-11
//         "?ptnr1_auth_seq_id", "?ptnr2_auth_seq_id",             // 12-13
//         "?ptnr1_label_seq_id", "?ptnr2_label_seq_id",           // 14-15
//         "?pdbx_ptnr1_PDB_ins_code", "?pdbx_ptnr2_PDB_ins_code", // 16-17
//         "?ptnr1_symmetry", "?ptnr2_symmetry",                   // 18-19
//         "?pdbx_dist_value", "?ccp4_link_id"})) {                // 20-21
//     Connection c;
//     c.name = row.str(kId);
//     copy_string(row, kLinkId, c.link_id);
//     c.type = connection_type_from_string(row.str(kConnTypeId));
//     if (row.has2(kSym1) && row.has2(kSym2)) {
//       c.asu = (row.str(kSym1) == row.str(kSym2) ? Asu::Same : Asu::Different);
//     }
//     copy_double(row, kDistValue, c.reported_distance);
//     for (int i = 0; i < 2; ++i) {
//       AtomAddress& a = (i == 0 ? c.partner1 : c.partner2);
//       if (row.has(kAuthAsymId+i) && row.has(kAuthSeqId+i)) {
//         a.chain_name = row.str(kAuthAsymId+i);
//         a.res_id = make_resid(row.str(kLabelCompId+i),
//                               row.str(kAuthSeqId+i), row.ptr_at(kInsCode+i));
//       } else if (row.has(kLabelAsymId+i) && row.has(kLabelSeqId+i)) {
//         set_part_of_address_from_label(a, st.first_model(),
//                                        row.str(kLabelAsymId+i),
//                                        row[kLabelSeqId+i]);
//         a.res_id.name = row.str(kLabelCompId+i);
//       } else {
//         fail("_struct_conn without either _auth_ or _label_ asym_id+seq_id");
//       }
//       a.atom_name = row.str(kLabelAtomId+i);
//       if (row.has2(kLabelAltId+i))
//         a.altloc = cif::as_char(row[kLabelAltId+i], '\0');
//     }
//     st.connections.emplace_back(c);
//   }
// }

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
