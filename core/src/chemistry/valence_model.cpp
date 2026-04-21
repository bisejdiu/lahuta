/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr std::string_view first = "besian", last = "sejdiu", host = "gmail.com";
 *   return std::string(first) + std::string(last) + "@" + std::string(host);
 * }();
 *
 */

#include <rdkit/GraphMol/MonomerInfo.h>

#include "chemistry/utils.hpp"
#include "chemistry/valence_model.hpp"
#include "residues/definitions.hpp"

namespace lahuta {

namespace {

enum class HydrogenScopeKind {
  Polymer,
  Water,
  Other
};

HydrogenScopeKind classify_scope_kind(const Residue &residue) {
  if (definitions::is_water(residue.name)) return HydrogenScopeKind::Water;
  if (definitions::is_polymer(residue.name)) return HydrogenScopeKind::Polymer;
  return HydrogenScopeKind::Other;
}

bool residue_has_explicit_hydrogen(const Residue &residue) {
  for (const auto *atom : residue.atoms) {
    if (atom && atom->getAtomicNum() == 1) return true;
  }
  return false;
}

} // namespace

bool ValenceModel::has_double_or_aromatic_bond(const RDKit::Atom &atom, const RDKit::ROMol &mol) const {
  for (const auto &bond : mol.atomBonds(&atom)) {
    auto *other_atom = bond->getOtherAtom(&atom);
    if (is_metal_atom(*other_atom)) continue;
    auto bond_type = bond->getBondType();
    if (bond_type == RDKit::Bond::DOUBLE || bond_type == RDKit::Bond::AROMATIC) {
      return true;
    }
  }
  return false;
}

bool ValenceModel::is_excluded_bond(const RDKit::Atom &atom_b, const RDKit::Atom &atom_c) const {
  unsigned int atomic_num_b = atom_b.getAtomicNum();
  unsigned int atomic_num_c = atom_c.getAtomicNum();
  return ((atomic_num_b == 15 || atomic_num_b == 16) && atomic_num_c == 8);
}

bool ValenceModel::is_conjugated_and_not_excluded(const RDKit::Atom &atom_b, const RDKit::ROMol &mol) const {
  for (const auto &bond_b : mol.atomBonds(&atom_b)) {
    auto *other_atom_b = bond_b->getOtherAtom(&atom_b);
    if (is_metal_atom(*other_atom_b)) continue;
    auto bond_type = bond_b->getBondType();
    if (bond_type == RDKit::Bond::DOUBLE || bond_type == RDKit::Bond::AROMATIC) {
      RDKit::Atom *atom_c = other_atom_b;
      if (!is_excluded_bond(atom_b, *atom_c)) {
        return true;
      }
    }
  }
  return false;
}

bool ValenceModel::is_conjugated(const RDKit::ROMol &mol, const RDKit::Atom &atom) const {
  unsigned int atomic_num = atom.getAtomicNum();
  bool is_hetero          = (atomic_num == 7 || atomic_num == 8);

  if (is_hetero && atom.getDegree() == 4) {
    return false;
  }

  if (has_double_or_aromatic_bond(atom, mol)) {
    return true;
  }

  // If atom is N or O, inspect its directly bonded neighbors
  if (is_hetero) {
    for (const auto &bond : mol.atomBonds(&atom)) {
      const RDKit::Atom *neighbor = bond->getOtherAtom(&atom);
      if (is_metal_atom(*neighbor)) continue;
      if (neighbor && is_conjugated_and_not_excluded(*neighbor, mol)) {
        return true;
      }
    }
  }

  return false;
}

int ValenceModel::get_element_count(const RDKit::ROMol &mol, const RDKit::Atom &atom, int element) const {
  int count = 0;
  for (const auto &bond : mol.atomBonds(&atom)) {
    auto other_atom = bond->getOtherAtom(&atom);
    if (is_metal_atom(*other_atom)) continue;
    if (other_atom->getAtomicNum() == element) {
      count++;
    }
  }
  return count;
}

bool ValenceModel::is_bound_to_sulfur_or_metal(const RDKit::ROMol &mol, RDKit::Atom *atom) const {
  for (const auto &bond : mol.atomBonds(atom)) {
    RDKit::Atom *other_atom       = bond->getOtherAtom(atom);
    unsigned int other_atomic_num = other_atom->getAtomicNum();

    auto name   = other_atom->getSymbol();
    gemmi::El e = gemmi::find_element(name.c_str());

    if (other_atomic_num == 16 || gemmi::is_metal(e)) {
      return true;
    }
  }
  return false;
}

bool ValenceModel::is_amidine_or_guanidine_nitrogen(int degree, int h_count, int valence) const {
  return (degree - h_count == 1 && valence - h_count == 2);
}

std::vector<bool> ValenceModel::compute_scope_has_explicit_h(const RDKit::RWMol &mol, const Residues &residues) const {
  std::vector<bool> atom_scope_has_explicit_h(mol.getNumAtoms(), false);

  const auto &residue_list = residues.get_residues();
  std::size_t i = 0;
  while (i < residue_list.size()) {
    const auto &seed = residue_list[i];
    const auto kind = classify_scope_kind(seed);

    std::size_t j = i + 1;
    switch (kind) {
      case HydrogenScopeKind::Polymer:
        while (j < residue_list.size()) {
          const auto &next = residue_list[j];
          if (classify_scope_kind(next) != HydrogenScopeKind::Polymer) break;
          if (next.chain_id != seed.chain_id) break;
          ++j;
        }
        break;
      case HydrogenScopeKind::Water:
        while (j < residue_list.size()) {
          const auto &next = residue_list[j];
          if (classify_scope_kind(next) != HydrogenScopeKind::Water) break;
          if (next.chain_id != seed.chain_id) break;
          ++j;
        }
        break;
      case HydrogenScopeKind::Other:
        break;
    }

    bool scope_has_explicit_h = false;
    for (std::size_t k = i; k < j; ++k) {
      if (residue_has_explicit_hydrogen(residue_list[k])) {
        scope_has_explicit_h = true;
        break;
      }
    }

    for (std::size_t k = i; k < j; ++k) {
      for (const auto *atom : residue_list[k].atoms) {
        if (!atom) continue;
        atom_scope_has_explicit_h[atom->getIdx()] = scope_has_explicit_h;
      }
    }

    i = j;
  }

  return atom_scope_has_explicit_h;
}

bool ValenceModel::has_neighbor_with_double_bonded_oxygen(const RDKit::ROMol &mol, RDKit::Atom *atom) const {
  for (const auto &bond_a : mol.atomBonds(atom)) {
    RDKit::Atom *other_atom_a = bond_a->getOtherAtom(atom);
    if (is_metal_atom(*other_atom_a)) continue;

    for (const auto &bond_b : mol.atomBonds(other_atom_a)) {
      RDKit::Atom *other_atom_b = bond_b->getOtherAtom(other_atom_a);
      if (is_metal_atom(*other_atom_b)) continue;

      if (other_atom_b != atom && other_atom_b->getAtomicNum() == 8 &&
          bond_b->getBondType() == RDKit::Bond::DOUBLE) {
        return true;
      }
    }
  }
  return false;
}

HybridizationType ValenceModel::assign_geometry(int total_coordination) const {
  switch (total_coordination) {
    case 0:
      return HybridizationType::S;
    case 1:
      return HybridizationType::S;
    case 2:
      return HybridizationType::SP;
    case 3:
      return HybridizationType::SP2;
    case 4:
      return HybridizationType::SP3;
    default:
      return HybridizationType::UNSPECIFIED;
  }
}

// FIX: TRP-NE1 is not evaluated as having 1 implicit H (but is by PropKa)
//      HIS-NE1 is not evaluated as having 1 implicit H (but is by PropKa)
//
// NOTE:
//      ASP-OD2 is evaluated as having 1 implicit H by PropKa (at 7.4) but only
//      those facing the core of the protein. This means H assignment by PropKa
//      is context dependent (more accurate).
void ValenceModel::molstar_valence_model(const RDKit::ROMol &mol, RDKit::Atom &atom, bool scope_has_explicit_h) {

  const int h_count             = get_element_count(mol, atom, 1);
  const unsigned int atomic_num = atom.getAtomicNum();
  int charge                    = atom.getFormalCharge();

  const bool assign_charge_flag = charge == 0;
  const bool assign_h_flag      = h_count == 0 && !scope_has_explicit_h;

  const int metal_bonds = static_cast<int>(count_bonds_to_metals(mol, atom));
  const int degree  = static_cast<int>(atom.getDegree()) - metal_bonds;
  const int valence = atom.getExplicitValence() - metal_bonds;

  const bool conjugated = is_conjugated(mol, atom);
  const bool multi_bond = (valence - degree > 0);

  int implicit_h = 0;
  auto geom      = HybridizationType::UNSPECIFIED;
  const auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom.getMonomerInfo());

  switch (atomic_num) {
    case 1:
      if (assign_charge_flag) {
        atom.setHybridization(HybridizationType::S);
        if (degree == 0) {
          // Hydrogen with no bonds has a positive charge
          charge = 1;
        } else if (degree == 1) {
          // Hydrogen with one bond is neutral
          charge = 0;
        }
      }
      // FIX: remove parentheses
      {
        atom.setFormalCharge(charge);
      }
      break;

    case 6:
      if (assign_charge_flag) {
        charge = 0; // typically neutral
        atom.setFormalCharge(charge);
      }
      if (assign_h_flag) {
        // normally forms 4 bonds: adjust implicit hydrogens
        implicit_h = std::max(0, 4 - valence - std::abs(charge));
      }
      {
        // Carbocation is planar, carbanion is tetrahedral
        int total_bonds = degree + implicit_h + std::max(0, -charge);
        atom.setHybridization(assign_geometry(total_bonds));
      }
      break;

    case 7:
      if (assign_charge_flag) {
        if (!assign_h_flag) {
          // Trust input H-count and calculate charge
          charge = valence - 3;
        } else if (conjugated && valence < 4) {
          // Conjugated nitrogen typically neutral unless special cases
          if (is_amidine_or_guanidine_nitrogen(degree, h_count, valence)) {
            charge = 1; // Example: Amidine or guanidine N
          } else {
            charge = 0; // Default to neutral
          }
        } else {
          // Check for special cases like sulfonamide N or metals
          if (is_bound_to_sulfur_or_metal(mol, &atom)) {
            charge = 0; // Nitrogen bound to sulfur or metal is neutral
          } else {
            charge = 1; // Otherwise assume positive charge
          }
        }
      }

      if (assign_h_flag) {
        // NH4+ -> 4, 1' amide -> 2, nitro N/N+ depiction -> 0
        implicit_h = std::max(0, 3 - valence + charge);
      }

      if (conjugated && !multi_bond) {
        // Conjugated nitrogen (amide, anilinic N) geometry treated differently
        // Assuming trigonal geometry for simplicity
        geom = assign_geometry(degree + implicit_h - charge);
      } else {
        // Normal nitrogen geometry (pyridine, amine, nitrile)
        geom = assign_geometry(degree + implicit_h + 1 - charge);
      }
      atom.setHybridization(geom);
      atom.setFormalCharge(charge);
      break;

    case 8:
      if (assign_charge_flag) {
        if (!assign_h_flag) {
          charge = valence - 2;
        }

        if (valence == 1 && has_neighbor_with_double_bonded_oxygen(mol, &atom)) {
          charge = -1;
        }
      }

      if (assign_h_flag) {
        implicit_h = std::max(0, 2 - valence + charge);
      }

      geom = assign_geometry(degree + implicit_h - charge + (conjugated && !multi_bond ? 1 : 2));

      atom.setHybridization(static_cast<HybridizationType>(geom));
      atom.setFormalCharge(charge);
      break;

    case 16:
      if (assign_charge_flag) {
        if (!assign_h_flag) {
          if (valence <= 3 && h_count == 0) {
            charge = valence - 2; // Deprotonated thiol
          } else {
            charge = 0;
          }
          atom.setFormalCharge(charge);
        }
      }
      if (assign_h_flag) {
        if (valence < 2) {
          implicit_h = std::max(0, 2 - valence + charge);
        }
      }
      if (valence <= 3) {
        // Thiol, thiolate, tioether -> tetrahedral
        geom = assign_geometry(degree + implicit_h - charge + 2);
        atom.setHybridization(geom);
      }

      break;

    // Halogens
    case 9:  // Fluorine (F)
    case 17: // Chlorine (Cl)
    case 35: // Bromine (Br)
    case 53: // Iodine (I)
    case 85: // Astatine (At)
      // Never implicitly protonate halides
      if (assign_charge_flag) {
        atom.setFormalCharge(valence - 1);
      }
      break;

    // Alkali metals
    case 3:  // Lithium (Li)
    case 11: // Sodium (Na)
    case 19: // Potassium (K)
    case 37: // Rubidium (Rb)
    case 55: // Cesium (Cs)
    case 87: // Francium (Fr)
      if (assign_charge_flag) {
        atom.setFormalCharge(1 - valence);
      }
      break;

    // Alkaline earth metals
    case 4:  // Beryllium (Be)
    case 12: // Magnesium (Mg)
    case 20: // Calcium (Ca)
    case 38: // Strontium (Sr)
    case 56: // Barium (Ba)
    case 88: // Radium (Ra)
      if (assign_charge_flag) {
        atom.setFormalCharge(2 - valence);
      }
      break;
    default:
      break;
  }

  // FIX: temporary, until we debug the issue
  if (assign_h_flag && res_info->getResidueName() == "TRP" && res_info->getName() == "NE1") {
    atom.setNumCompImplicitHs(1);
    return;
  } else if (assign_h_flag && res_info->getResidueName() == "HIS" && res_info->getName() == "NE2") {
    atom.setNumCompImplicitHs(1);
    return;
  }

  atom.setNumCompImplicitHs(implicit_h);
}

} // namespace lahuta
