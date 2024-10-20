#include "contacts/charges.hpp"
#include "contacts/heuristics.hpp"
#include "contacts/utils.hpp"
#include <unordered_map>

namespace lahuta {

Feature create_feature(AtomType type, FeatureGroup group, const std::vector<const RDKit::Atom *> &members) {
  Feature feature;
  feature.type = type;
  feature.group = group;
  feature.members = members;
  return feature;
}

auto identify_feature_groups(const RDKit::RWMol &mol) {
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

auto identify_negative_feature_groups(const RDKit::RWMol &mol) {
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

std::vector<Feature> add_positive_charges(const RDKit::RWMol &mol, ResMap &res_map) {

  std::vector<Feature> features;
  /*IncrementalVector features;*/
  std::optional<std::unordered_map<const RDKit::Atom *, FeatureGroup>> groups;
  std::unordered_set<const RDKit::Atom *> added_atoms;

  for (const auto &residue_pair : res_map) {
    auto &atoms_in_residue = residue_pair.second;
    auto *res_info =
        static_cast<const RDKit::AtomPDBResidueInfo *>(atoms_in_residue.front()->getMonomerInfo());
    std::string res_name = res_info->getResidueName();

    // Handle positively charged residues (ARG, HIS, LYS)
    if (PositivelyChargedResidues.count(res_name)) {
      auto feature = create_feature(AtomType::POS_IONISABLE, FeatureGroup::None, {});

      for (const auto *atom : atoms_in_residue) {
        auto *atom_res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
        std::string atom_name = atom_res_info->getName();

        if (atom->getAtomicNum() == 7 && ProteinBackboneAtoms.count(atom_name) == 0) {
          feature.members.push_back(atom);
          added_atoms.insert(atom);
        }
      }

      if (!feature.members.empty()) {
        features.push_back(feature);
      }
    } else if (PolymerNames.count(res_name) == 0) {
      // Handle non-polymer residues
      if (!groups.has_value()) {
        groups = identify_feature_groups(mol);
      }
      for (const auto *atom : atoms_in_residue) {
        auto it = groups->find(atom);
        if (it != groups->end()) {
          auto nitrogens = bonded_atoms(mol, atom, 7);
          if (!nitrogens.empty()) {
            auto feature = create_feature(AtomType::POS_IONISABLE, it->second, nitrogens);
            features.push_back(feature);
            added_atoms.insert(nitrogens.begin(), nitrogens.end());
          }
        } else {
          // Add remaining positively charged atoms not already added
          if (atom->getFormalCharge() > 0 && added_atoms.count(atom) == 0) {
            auto feature = create_feature(AtomType::POS_IONISABLE, FeatureGroup::None, {atom});
            features.push_back(feature);
            added_atoms.insert(atom);
          }
        }
      }
    }
  }

  return std::move(features);
}

std::vector<Feature> add_negative_charges(const RDKit::RWMol &mol, ResMap &res_map) {

  std::optional<std::unordered_map<const RDKit::Atom *, FeatureGroup>> groups;
  std::vector<Feature> features;
  /*IncrementalVector features;*/
  std::unordered_set<const RDKit::Atom *> added_atoms;

  for (const auto &residue_pair : res_map) {
    const std::vector<const RDKit::Atom *> &atoms_in_residue = residue_pair.second;
    auto *res_info =
        static_cast<const RDKit::AtomPDBResidueInfo *>(atoms_in_residue.front()->getMonomerInfo());
    std::string res_name = res_info->getResidueName();

    if (NegativelyChargedResidues.count(res_name)) {
      // Handle negatively charged residues (GLU, ASP)
      Feature feature;
      feature.type = AtomType::NEG_IONISABLE;
      feature.group = FeatureGroup::None;

      for (const auto *atom : atoms_in_residue) {
        auto *atom_res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
        std::string atom_name = atom_res_info->getName();

        if (atom->getAtomicNum() == 8 && ProteinBackboneAtoms.count(atom_name) == 0) {
          feature.members.push_back(atom);
          added_atoms.insert(atom);
        }
      }

      if (!feature.members.empty()) {
        features.push_back(feature);
      }
    } else if (BaseNames.count(res_name)) {
      // Handle nucleic acid bases (DNA/RNA)
      for (const auto *atom : atoms_in_residue) {
        if (is_phosphate(mol, *atom)) {
          auto oxygens = bonded_atoms(mol, atom, 8);
          if (!oxygens.empty()) {
            auto feature = create_feature(AtomType::NEG_IONISABLE, FeatureGroup::Phosphate, oxygens);
            features.push_back(feature);
            added_atoms.insert(oxygens.begin(), oxygens.end());
          }
        }
      }
    } else if (PolymerNames.count(res_name) == 0) {
      // Handle non-polymer residues
      if (!groups.has_value()) {
        groups = identify_negative_feature_groups(mol);
      }
      for (const auto *atom : atoms_in_residue) {
        auto it = groups->find(atom);
        if (it != groups->end()) {
          auto oxygens = bonded_atoms(mol, atom, 8);
          if (!oxygens.empty()) {
            auto feature = create_feature(AtomType::NEG_IONISABLE, it->second, oxygens);
            features.push_back(feature);
            added_atoms.insert(oxygens.begin(), oxygens.end());
          }
        } else {
          // Collect non-backbone nitrogen atoms
          auto *atom_res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
          std::string atom_name = atom_res_info->getName();

          if (atom->getAtomicNum() == 7 && ProteinBackboneAtoms.count(atom_name) == 0) {
            auto feature = create_feature(AtomType::NEG_IONISABLE, FeatureGroup::None, {atom});
            features.push_back(feature);
            added_atoms.insert(atom);
          }
        }

        // Add remaining negatively charged atoms not already added
        if (atom->getFormalCharge() < 0 && added_atoms.count(atom) == 0) {
          auto feature = create_feature(AtomType::NEG_IONISABLE, FeatureGroup::None, {atom});
          features.push_back(feature);
          added_atoms.insert(atom);
        }
      }
    }
  }

  return std::move(features);
}

} // namespace lahuta
