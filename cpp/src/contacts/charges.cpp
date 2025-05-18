#include "contacts/charges.hpp"
#include "contacts/heuristics.hpp"
#include "contacts/utils.hpp"
#include "definitions.hpp"
#include "residues.hpp"
#include <elements.hpp>
#include <unordered_map>

// clang-format off
namespace lahuta {

auto identify_positive_charge_groups(const RDKit::RWMol &mol) {
  std::unordered_map<const RDKit::Atom *, FeatureGroup> groups;

  for (const auto *atom : mol.atoms()) {
    FeatureGroup group = FeatureGroup::None;

    if      (is_C_in_guanidine  (mol, *atom)) { group = FeatureGroup::Guanidine; }
    else if (is_C_in_acetamidine(mol, *atom)) { group = FeatureGroup::Acetamidine; }

    if (group != FeatureGroup::None) {
      groups[atom] = group;
    }
  }

  return groups;
}

auto identify_negative_charge_groups(const RDKit::RWMol &mol) {
  std::unordered_map<const RDKit::Atom *, FeatureGroup> groups;

  for (const auto *atom : mol.atoms()) {
    FeatureGroup group = FeatureGroup::None;

    if      (is_S_in_sulfonic_acid(mol, *atom)) { group = FeatureGroup::SulfonicAcid; }
    else if (is_P_in_phosphate    (mol, *atom)) { group = FeatureGroup::Phosphate; }
    else if (is_S_in_sulfate      (mol, *atom)) { group = FeatureGroup::Sulfate; }
    else if (is_C_in_carboxylate  (mol, *atom)) { group = FeatureGroup::Carboxylate; }

    if (group != FeatureGroup::None) {
      groups[atom] = group;
    }
  }

  return groups;
}

GroupEntityCollection add_positive_charges(const RDKit::RWMol &mol, const Residues &residues) {

  std::optional<std::unordered_map<const RDKit::Atom *, FeatureGroup>> groups;
  std::unordered_set<const RDKit::Atom *> added_atoms;

  GroupEntityCollection features;
  for (const auto &residue : residues) {

    // Handle positively charged residues (ARG, HIS, LYS)
    if (definitions::is_positive_charge(residue.name)) {
      std::vector<const RDKit::Atom *> members;
      for (const auto *atom : residue.atoms) {
        auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
        std::string atom_name = res_info->getName();

        if (atom->getAtomicNum() == Element::N && definitions::ProteinBackboneAtoms.count(atom_name) == 0) {
          members.push_back(atom);
          added_atoms.insert(atom);
        }
      }

      if (!members.empty()) {
        features.add_data(AtomType::POS_IONISABLE, FeatureGroup::None, members);
      }
    } else if (!definitions::is_polymer(residue.name)) {
      // Handle non-polymer residues
      if (!groups.has_value()) {
        groups = identify_positive_charge_groups(mol);
      }

      for (const auto *atom : residue.atoms) {
        auto it = groups->find(atom);
        if (it != groups->end()) {
          auto nitrogens = bonded_atoms(mol, atom, Element::N);
          if (nitrogens.empty()) continue;

          features.add_data(AtomType::POS_IONISABLE, it->second, nitrogens);
          added_atoms.insert(nitrogens.begin(), nitrogens.end());
        } else {
          // Add remaining positively charged atoms not already added
          if (atom->getFormalCharge() > 0 && added_atoms.count(atom) == 0) {
            features.add_data(AtomType::POS_IONISABLE, FeatureGroup::None, {atom});
            added_atoms.insert(atom);
          }
        }
      }
    }
  }

  return std::move(features);
}

GroupEntityCollection add_negative_charges(const RDKit::RWMol &mol, const Residues &residues) {

  std::optional<std::unordered_map<const RDKit::Atom *, FeatureGroup>> groups;
  std::unordered_set<const RDKit::Atom *> added_atoms;

  GroupEntityCollection features;
  for (const auto &residue : residues) {
    if (definitions::is_negative_charge(residue.name)) {

      // Handle negatively charged residues (GLU, ASP)
      std::vector<const RDKit::Atom *> members;
      for (const auto *atom : residue.atoms) {
        auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
        std::string atom_name = res_info->getName();

        if (atom->getAtomicNum() == Element::O && definitions::ProteinBackboneAtoms.count(atom_name) == 0) {
          members.push_back(atom);
          added_atoms.insert(atom);
        }
      }

      if (members.empty()) continue;
      features.add_data(AtomType::NEG_IONISABLE, FeatureGroup::None, members);
    } else if (definitions::is_base(residue.name)) {
      // Handle nucleic acid bases (DNA/RNA)
      for (const auto *atom : residue.atoms) {
        if (is_P_in_phosphate(mol, *atom)) {
          auto oxygens = bonded_atoms(mol, atom, Element::O);
          if (oxygens.empty()) continue;

          features.add_data(AtomType::NEG_IONISABLE, FeatureGroup::Phosphate, oxygens);
          added_atoms.insert(oxygens.begin(), oxygens.end());
        }
      }
    } else if (!definitions::is_polymer(residue.name)) {
      // Handle non-polymer residues
      if (!groups.has_value()) {
        groups = identify_negative_charge_groups(mol);
      }
      for (const auto *atom : residue.atoms) {
        auto it = groups->find(atom);
        if (it != groups->end()) {
          auto oxygens = bonded_atoms(mol, atom, Element::O);
          if (oxygens.empty()) continue;

          features.add_data(AtomType::NEG_IONISABLE, it->second, oxygens);
          added_atoms.insert(oxygens.begin(), oxygens.end());

        } else {
          // Collect non-backbone nitrogen atoms
          auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
          std::string atom_name = res_info->getName();

          if (atom->getAtomicNum() == Element::N && definitions::ProteinBackboneAtoms.count(atom_name) == 0) {
            features.add_data(AtomType::NEG_IONISABLE, FeatureGroup::None, {atom});
            added_atoms.insert(atom);
          }
        }

        // Add remaining negatively charged atoms not already added
        if (atom->getFormalCharge() < 0 && added_atoms.count(atom) == 0) {
          features.add_data(AtomType::NEG_IONISABLE, FeatureGroup::None, {atom});
          added_atoms.insert(atom);
        }
      }
    }
  }

  return std::move(features);
}

} // namespace lahuta
