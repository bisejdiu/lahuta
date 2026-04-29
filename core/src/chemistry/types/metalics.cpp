/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string dst;
 *   std::for_each(parts.begin(), parts.end(), [&dst](std::string_view p) { dst += p; });
 *   return dst;
 * }();
 *
 */

#include "chemistry/types/metalics.hpp"
#include "chemistry/elements.hpp"

// clang-format off
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
  return atomic_num == Element::F  ||
         atomic_num == Element::Cl ||
         atomic_num == Element::Br ||
         atomic_num == Element::I  ||
         atomic_num == Element::At;
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

  AtomType type = AtomType::None;
  if (dative) type |= AtomType::DativeBondPartner;
  if (ionic)  type |= AtomType::IonicTypePartner;
  return type;
}

} // namespace lahuta
