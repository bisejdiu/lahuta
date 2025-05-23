#include "contacts/metals.hpp"

namespace lahuta {

bool is_transition_metal(int atomic_num) {
  return (atomic_num >= 21 && atomic_num <= 29) || // Ti to Cu
         (atomic_num >= 39 && atomic_num <= 47) || // Y to Ag
         (atomic_num >= 72 && atomic_num <= 79) || // Hf to Au
         (atomic_num >= 104 && atomic_num <= 108); // Rf to Hs
}

AtomType add_metal(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  const int atomic_num = atom.getAtomicNum();
  AtomType type = AtomType::None;

  if (std::find(IonicTypeMetals.begin(), IonicTypeMetals.end(), atomic_num) != IonicTypeMetals.end()) {
    type = AtomType::IonicTypeMetal;
  } else if (is_transition_metal(atomic_num) || atomic_num == 30 || atomic_num == 48) {
    type = AtomType::TransitionMetal;
  }

  return type;
}

bool is_protein_sidechain(const std::string &atomname) {
  return definitions::ProteinBackboneAtoms.find(atomname) == definitions::ProteinBackboneAtoms.end();
}

bool is_protein_backbone(const std::string &atomname) {
  return definitions::ProteinBackboneAtoms.find(atomname) != definitions::ProteinBackboneAtoms.end();
}

bool is_nucleic_backbone(const std::string &atomname) {
  return definitions::NucleicBackboneAtoms.find(atomname) != definitions::NucleicBackboneAtoms.end();
}

bool is_halogen(int atomic_num) {
  // Halogens: F(9), Cl(17), Br(35), I(53), At(85)
  return atomic_num == 9 || atomic_num == 17 || atomic_num == 35 || atomic_num == 53 || atomic_num == 85;
}

AtomType add_metal_binding(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  int atomic_num = atom.getAtomicNum();

  std::string resname;
  std::string atomname;
  const auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom.getMonomerInfo());

  if (!res_info) {
    return AtomType::None;
  }

  resname = res_info->getResidueName();
  atomname = res_info->getName();

  bool dative = false, ionic = false;

  bool is_standard_aminoacid = definitions::is_protein_extended(resname);
  bool is_standard_base = definitions::is_base(resname);

  if (!is_standard_aminoacid && !is_standard_base) {
    if (is_halogen(atomic_num) || atomic_num == 8 || atomic_num == 16) {
      dative = true, ionic = true;
    } else if (atomic_num == 7) {
      dative = true;
    }
  } else if (is_standard_aminoacid) {
    // Main chain oxygen atom or oxygen, nitrogen, and sulfur from specific amino acids
    if (atomic_num == 8) {
      if (is_O_metal_binding(resname) && is_protein_sidechain(atomname)) {
        dative = true, ionic = true;
      } else if (is_protein_backbone(atomname)) {
        dative = true, ionic = true;
      }
    } else if (atomic_num == 16 && (resname == "CYS" || resname == "MET")) {
      dative = true, ionic = true;
    } else if (atomic_num == 7) {
      if (resname == "HIS" && is_protein_sidechain(atomname)) {
        dative = true;
      }
    }
  } else if (is_standard_base) {
    if (atomic_num == 8 && is_nucleic_backbone(atomname)) {
      dative = true, ionic = true;
    } else if (atomname == "N3" || atomname == "N4" || atomname == "N7") {
      dative = true;
    } else if (atomname == "O2" || atomname == "O4" || atomname == "O6") {
      dative = true, ionic = true;
    }
  }

  if (dative) {
    return AtomType::DativeBondPartner;
  }
  if (ionic) {
    return AtomType::IonicTypePartner;
  }
  return AtomType::None;
}

} // namespace lahuta
