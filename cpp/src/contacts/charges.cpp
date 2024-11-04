#include "contacts/charges.hpp"
#include "contacts/heuristics.hpp"
#include "contacts/utils.hpp"
#include "definitions.hpp"
#include "residues.hpp"
#include <unordered_map>

namespace lahuta {

auto identify_positive_charge_groups(const RDKit::RWMol &mol) {
  std::unordered_map<const RDKit::Atom *, FeatureGroup> identified_atoms;

  for (const auto *atom : mol.atoms()) {
    FeatureGroup group = FeatureGroup::None;

    if (is_guanidine(mol, *atom)) {
      group = FeatureGroup::Guanidine;
    } else if (is_acetamidine(mol, *atom)) {
      group = FeatureGroup::Acetamidine;
    }

    if (group != FeatureGroup::None) {
      identified_atoms[atom] = group;
    }
  }

  return identified_atoms;
}

auto identify_negative_charge_groups(const RDKit::RWMol &mol) {
  std::unordered_map<const RDKit::Atom *, FeatureGroup> identified_atoms;

  for (const auto *atom : mol.atoms()) {
    FeatureGroup group = FeatureGroup::None;

    if (is_sulfonic_acid(mol, *atom)) {
      group = FeatureGroup::SulfonicAcid;
    } else if (is_phosphate(mol, *atom)) {
      group = FeatureGroup::Phosphate;
    } else if (is_sulfate(mol, *atom)) {
      group = FeatureGroup::Sulfate;
    } else if (is_carboxylate(mol, *atom)) {
      group = FeatureGroup::Carboxylate;
    }

    if (group != FeatureGroup::None) {
      identified_atoms[atom] = group;
    }
  }

  return identified_atoms;
}

GroupEntityCollection add_positive_charges(const RDKit::RWMol &mol, const Residues &residues) {

  std::optional<std::unordered_map<const RDKit::Atom *, FeatureGroup>> groups;
  std::unordered_set<const RDKit::Atom *> added_atoms;

  GroupEntityCollection features;
  for (const auto &residue : residues) {

    // Handle positively charged residues (ARG, HIS, LYS)
    if (definitions::PositivelyChargedResidues.count(residue.name)) {

      std::vector<const RDKit::Atom *> members;
      for (const auto *atom : residue.atoms) {
        auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
        std::string atom_name = res_info->getName();

        if (atom->getAtomicNum() == 7 && definitions::ProteinBackboneAtoms.count(atom_name) == 0) {
          members.push_back(atom);
          added_atoms.insert(atom);
        }
      }

      if (!members.empty()) {
        features.add_data(AtomType::POS_IONISABLE, FeatureGroup::None, members);
      }
    } else if (definitions::PolymerNames.count(residue.name) == 0) {
      // Handle non-polymer residues
      if (!groups.has_value()) {
        groups = identify_positive_charge_groups(mol);
      }

      for (const auto *atom : residue.atoms) {
        auto it = groups->find(atom);
        if (it != groups->end()) {
          auto nitrogens = bonded_atoms(mol, atom, 7);
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
    if (definitions::NegativelyChargedResidues.count(residue.name)) {

      // Handle negatively charged residues (GLU, ASP)
      std::vector<const RDKit::Atom *> members;
      for (const auto *atom : residue.atoms) {
        auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
        std::string atom_name = res_info->getName();

        if (atom->getAtomicNum() == 8 && definitions::ProteinBackboneAtoms.count(atom_name) == 0) {
          members.push_back(atom);
          added_atoms.insert(atom);
        }
      }

      if (members.empty()) continue;
      features.add_data(AtomType::NEG_IONISABLE, FeatureGroup::None, members);
    } else if (definitions::BaseNames.count(residue.name)) {
      // Handle nucleic acid bases (DNA/RNA)
      for (const auto *atom : residue.atoms) {
        if (is_phosphate(mol, *atom)) {
          auto oxygens = bonded_atoms(mol, atom, 8);
          if (oxygens.empty()) continue;

          features.add_data(AtomType::NEG_IONISABLE, FeatureGroup::Phosphate, oxygens);
          added_atoms.insert(oxygens.begin(), oxygens.end());
        }
      }
    } else if (definitions::PolymerNames.count(residue.name) == 0) {
      // Handle non-polymer residues
      if (!groups.has_value()) {
        groups = identify_negative_charge_groups(mol);
      }
      for (const auto *atom : residue.atoms) {
        auto it = groups->find(atom);
        if (it != groups->end()) {
          auto oxygens = bonded_atoms(mol, atom, 8);
          if (oxygens.empty()) continue;

          features.add_data(AtomType::NEG_IONISABLE, it->second, oxygens);
          added_atoms.insert(oxygens.begin(), oxygens.end());

        } else {
          // Collect non-backbone nitrogen atoms
          auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
          std::string atom_name = res_info->getName();

          if (atom->getAtomicNum() == 7 && definitions::ProteinBackboneAtoms.count(atom_name) == 0) {
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
