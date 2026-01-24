#include <cctype>
#include <unordered_set>

#include <GraphMol/MonomerInfo.h>

#include "chemistry/elements.hpp"
#include "chemistry/types/hydrophobic.hpp"
#include "residues/definitions.hpp"
#include "typing/getcontacts/atom_typing.hpp"
#include "typing/types.hpp"

// clang-format off
namespace lahuta::typing::getcontacts {
namespace {

using namespace definitions;

inline bool has_bonded_element(const RDKit::RWMol& mol, const RDKit::Atom& atom, unsigned element) {
  for (const auto& bond : mol.atomBonds(&atom)) {
    if (bond->getOtherAtom(&atom)->getAtomicNum() == element) return true;
  }
  return false;
}

bool is_histidine_residue(const std::string& res) {
  return definitions::is_histidine(res);
}

bool is_positive_protein_site(const std::string& res, const std::string& name) {
  if (res == "LYS" && name == "NZ") return true;
  if (res == "ARG" && (name == "NH1" || name == "NH2")) return true;
  if (is_histidine_residue(res) && (name == "ND1" || name == "NE2")) return true;
  return false;
}

bool is_negative_protein_site(const std::string& res, const std::string& name) {
  if (res == "ASP" && (name == "OD1" || name == "OD2")) return true;
  if (res == "GLU" && (name == "OE1" || name == "OE2")) return true;

  static const std::unordered_set<std::string> nuc_residues = {"A","C","G","U","DA","DC","DG","DT"};
  if (nuc_residues.count(res) && (name == "OP1" || name == "OP2" || name == "O1P" || name == "O2P")) {
    return true;
  }
  return false;
}

bool is_carboxylate_oxygen(const RDKit::RWMol& mol, const RDKit::Atom& atom) {
  if (atom.getAtomicNum() != Element::O) return false;
  for (const auto& bond : mol.atomBonds(&atom)) {
    const auto* neighbor = bond->getOtherAtom(&atom);
    if (neighbor->getAtomicNum() != Element::C) continue;

    unsigned oxygen_count    = 0;
    unsigned carbon_count    = 0;
    unsigned neighbour_total = 0;
    for (const auto& cbond : mol.atomBonds(neighbor)) {
      const auto* other = cbond->getOtherAtom(neighbor);
      if (other == &atom) continue;
      ++neighbour_total;
      const auto z = other->getAtomicNum();
      if      (z == Element::O) ++oxygen_count;
      else if (z == Element::C) ++carbon_count;
    }

    if (oxygen_count >= 1 && (oxygen_count + carbon_count) >= 2 && neighbour_total >= 2) {
      return true;
    }
  }
  return false;
}

bool is_phosphate_oxygen(const RDKit::RWMol& mol, const RDKit::Atom& atom) {
  if (atom.getAtomicNum() != Element::O) return false;
  for (const auto& bond : mol.atomBonds(&atom)) {
    if (bond->getOtherAtom(&atom)->getAtomicNum() == Element::P) return true;
  }
  return false;
}

bool is_hbond_donor(const RDKit::RWMol& mol, const RDKit::Atom& atom, const std::string& res, const std::string& name) {
  const int atomic_num = atom.getAtomicNum();

  //
  // Backbone nitrogen is a potential donor (like VMD's measure hbonds)
  // This allows N-N backbone pairs to be detected as H-bonds
  // BUT exclude proline - it's a secondary amine with no H on backbone N
  //
  if (atomic_num == Element::N && name == "N" && res != "PRO") {
    return true;
  }

  if (atomic_num == Element::N || atomic_num == Element::O || atomic_num == Element::S) {
    // Check for explicit H atoms bonded OR implicit H count > 0
    // This allows donor detection to work with or without explicit hydrogens
    if (has_bonded_element(mol, atom, static_cast<unsigned int>(Element::H))) return true;
    if (atom.getTotalNumHs() > 0) return true;
  }

  if (is_histidine_residue(res) && (name == "ND1" || name == "NE2")) return true;
  return false;
}

bool is_hbond_acceptor(const RDKit::RWMol& mol, const RDKit::Atom& atom, const std::string& res, const std::string& name) {
  const int atomic_num = atom.getAtomicNum();

  if (atomic_num == Element::O) return true;

  if (atomic_num == Element::N) {
    //
    // Backbone nitrogen is always a acceptor (like VMD's measure hbonds)
    // This allows N-N backbone pairs to be detected as H-bonds. Proline CAN be an acceptor
    //
    if (name == "N") return true;

    if (is_histidine_residue(res) && (name == "ND1" || name == "NE2")) return true;
    if (atom.getFormalCharge() > 0) return false;

    const auto hybrid = atom.getHybridization();
    const unsigned neighbors = atom.getDegree() + atom.getNumImplicitHs();

    if ((hybrid == RDKit::Atom::HybridizationType::SP3 && neighbors < 4) ||
        (hybrid == RDKit::Atom::HybridizationType::SP2 && neighbors < 3) ||
        (hybrid == RDKit::Atom::HybridizationType::SP  && neighbors < 2)) {
      return true;
    }
  }

  if (atomic_num == Element::S) {
    if (res == "CYS" || res == "MET" || atom.getFormalCharge() < 0) return true;
  }

  if (atomic_num == Element::F) return true;
  return false;
}

bool is_ligand(const std::string& res) {
  if (res.empty()) return true;
  if (definitions::is_polymer(res))          return false;
  if (definitions::is_base(res))             return false;
  if (definitions::is_standard_protein(res)) return false;
  return true;
}

bool is_positive_ligand_site(const RDKit::Atom& atom) {
  static const std::unordered_set<unsigned> metal_cations = {
    static_cast<unsigned>(Element::Mg), static_cast<unsigned>(Element::Mn),
    static_cast<unsigned>(Element::Rh), static_cast<unsigned>(Element::Zn),
    static_cast<unsigned>(Element::Fe), static_cast<unsigned>(Element::Bi),
    static_cast<unsigned>(Element::As), static_cast<unsigned>(Element::Ag)
  };
  if (metal_cations.count(atom.getAtomicNum())) return true;
  if (atom.getFormalCharge() > 0) return true;
  return false;
}

bool is_negative_ligand_site(const RDKit::RWMol& mol, const RDKit::Atom& atom) {
  if (atom.getFormalCharge() < 0)       return true;
  if (is_carboxylate_oxygen(mol, atom)) return true;
  if (is_phosphate_oxygen(mol, atom))   return true;
  return false;
}

bool is_hydrophobic_residue(const std::string& res) {
  // Match getcontacts: only ALA, PHE, GLY, ILE, LEU, PRO, VAL, TRP, plus ligands, which are handled separately
  static const std::unordered_set<std::string> hydrophobic_residues = {"ALA", "PHE", "GLY", "ILE", "LEU", "PRO", "VAL", "TRP"};
  return hydrophobic_residues.count(res) > 0;
}

} // namespace

AtomType classify_atom(const RDKit::RWMol& mol, const RDKit::Atom& atom) {
  AtomType type = AtomType::None;

  //
  // Calling updatePropertyCache from the calling code does not seem to update
  // the code. I have no idea why. - Besian, October 2025
  //
  auto& atom_mut = const_cast<RDKit::Atom&>(atom);
  atom_mut.calcExplicitValence(false);
  atom_mut.calcImplicitValence(false);

  auto res_info = static_cast<const RDKit::AtomPDBResidueInfo*>(atom.getMonomerInfo());
  if (!res_info) return type;
  const auto res    = res_info->getResidueName();
  const auto name   = res_info->getName();
  const bool ligand = is_ligand(res);

  if (is_hbond_donor   (mol, atom, res, name)) type |= AtomType::HbondDonor;
  if (is_hbond_acceptor(mol, atom, res, name)) type |= AtomType::HbondAcceptor;

  if (!ligand && is_positive_protein_site(res, name)) type |= AtomType::PositiveCharge;
  if ( ligand && is_positive_ligand_site(atom))       type |= AtomType::PositiveCharge;

  if (!ligand && is_negative_protein_site(res, name)) type |= AtomType::NegativeCharge;
  if ( ligand && is_negative_ligand_site(mol, atom))  type |= AtomType::NegativeCharge;

  // Hydrophobic carbons / fluorine akin to getcontacts (carbon with only C/H neighbors).
  // Only for specific hydrophobic residues or ligands
  if (add_hydrophobic_atom(mol, atom) == AtomType::Hydrophobic) {
    if (ligand || is_hydrophobic_residue(res)) {
      type |= AtomType::Hydrophobic;
    }
  }

  if (atom.getIsAromatic()) type |= AtomType::Aromatic;

  return type;
}

} // namespace lahuta::typing::getcontacts
