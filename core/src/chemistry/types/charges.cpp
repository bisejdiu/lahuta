#include "chemistry/types/charges.hpp"
#include "chemistry/heuristics.hpp"
#include "chemistry/utils.hpp"
#include "definitions.hpp"
#include "entities/records.hpp"
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

std::vector<GroupRec> add_positive_charges(const RDKit::RWMol &mol, const Residues &residues) {
  std::optional<std::unordered_map<const RDKit::Atom *, FeatureGroup>> groups;
  std::unordered_set<const RDKit::Atom *> added_atoms;

  std::vector<GroupRec> features;
  for (const auto &residue : residues) {
    // Handle positively charged residues (ARG, HIS, LYS)
    if (definitions::is_positive_charge(residue.name)) {
      std::vector<std::reference_wrapper<const RDKit::Atom>> member_atoms;
      for (const auto *atom : residue.atoms) {
        auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
        std::string atom_name = res_info->getName();

        if (atom->getAtomicNum() == Element::N && definitions::ProteinBackboneAtoms.count(atom_name) == 0) {
          member_atoms.emplace_back(*atom);
          added_atoms.insert(atom);
        }
      }

      if (!member_atoms.empty()) {
        features.push_back(GroupRec{
          /*.a_type =*/ AtomType::PositiveCharge,
          /*.type   =*/ FeatureGroup::None,
          /*.atoms  =*/ std::move(member_atoms)
        });
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

          std::vector<std::reference_wrapper<const RDKit::Atom>> nitrogen_atoms;
          nitrogen_atoms.reserve(nitrogens.size());
          for (const auto* n_atom : nitrogens) {
            nitrogen_atoms.emplace_back(*n_atom);
            added_atoms.insert(n_atom);
          }

          features.push_back(GroupRec{
            /*.a_type =*/ AtomType::PositiveCharge,
            /*.type   =*/ it->second,
            /*.atoms  =*/ std::move(nitrogen_atoms)
          });
        } else {
          // Add remaining positively charged atoms not already added
          if (atom->getFormalCharge() > 0 && added_atoms.count(atom) == 0) {
            std::vector<std::reference_wrapper<const RDKit::Atom>> atom_refs = {std::cref(*atom)};
            features.push_back(GroupRec{
              /*.a_type =*/ AtomType::PositiveCharge,
              /*.type   =*/ FeatureGroup::None,
              /*.atoms  =*/ std::move(atom_refs)
            });
            added_atoms.insert(atom);
          }
        }
      }
    }
  }

  return features;
}

std::vector<GroupRec> add_negative_charges(const RDKit::RWMol &mol, const Residues &residues) {
  std::optional<std::unordered_map<const RDKit::Atom *, FeatureGroup>> groups;
  std::unordered_set<const RDKit::Atom *> added_atoms;

  std::vector<GroupRec> features;
  for (const auto &residue : residues) {
    if (definitions::is_negative_charge(residue.name)) {
      // Handle negatively charged residues (GLU, ASP)
      std::vector<std::reference_wrapper<const RDKit::Atom>> member_atoms;
      for (const auto *atom : residue.atoms) {
        auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
        std::string atom_name = res_info->getName();

        if (atom->getAtomicNum() == Element::O && definitions::ProteinBackboneAtoms.count(atom_name) == 0) {
          member_atoms.emplace_back(*atom);
          added_atoms.insert(atom);
        }
      }

      if (member_atoms.empty()) continue;

      features.push_back(GroupRec{
        /*.a_type =*/ AtomType::NegativeCharge,
        /*.type   =*/ FeatureGroup::None,
        /*.atoms  =*/ std::move(member_atoms)
      });
    } else if (definitions::is_base(residue.name)) {
      // Handle nucleic acid bases (DNA/RNA)
      for (const auto *atom : residue.atoms) {
        if (is_P_in_phosphate(mol, *atom)) {
          auto oxygens = bonded_atoms(mol, atom, Element::O);
          if (oxygens.empty()) continue;

          std::vector<std::reference_wrapper<const RDKit::Atom>> oxygen_atoms;
          oxygen_atoms.reserve(oxygens.size());
          for (const auto* o_atom : oxygens) {
            oxygen_atoms.emplace_back(*o_atom);
            added_atoms.insert(o_atom);
          }

          features.push_back(GroupRec{
            /*.a_type =*/ AtomType::NegativeCharge,
            /*.type   =*/ FeatureGroup::Phosphate,
            /*.atoms  =*/ std::move(oxygen_atoms)
          });
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

          std::vector<std::reference_wrapper<const RDKit::Atom>> oxygen_atoms;
          oxygen_atoms.reserve(oxygens.size());
          for (const auto* o_atom : oxygens) {
            oxygen_atoms.emplace_back(*o_atom);
            added_atoms.insert(o_atom);
          }

          features.push_back(GroupRec{
            /*.a_type =*/ AtomType::NegativeCharge,
            /*.type   =*/ it->second,
            /*.atoms  =*/ std::move(oxygen_atoms)
          });
        } else {
          // Collect non-backbone nitrogen atoms
          auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
          std::string atom_name = res_info->getName();

          if (atom->getAtomicNum() == Element::N && definitions::ProteinBackboneAtoms.count(atom_name) == 0) {
            std::vector<std::reference_wrapper<const RDKit::Atom>> atom_refs = {std::cref(*atom)};
            features.push_back(GroupRec{
              /*.a_type =*/ AtomType::NegativeCharge,
              /*.type   =*/ FeatureGroup::None,
              /*.atoms  =*/ std::move(atom_refs)
            });
            added_atoms.insert(atom);
          }
        }

        // Add remaining negatively charged atoms not already added
        if (atom->getFormalCharge() < 0 && added_atoms.count(atom) == 0) {
          std::vector<std::reference_wrapper<const RDKit::Atom>> atom_refs = {std::cref(*atom)};
          features.push_back(GroupRec{
            /*.a_type =*/ AtomType::NegativeCharge,
            /*.type   =*/ FeatureGroup::None,
            /*.atoms  =*/ std::move(atom_refs)
          });
          added_atoms.insert(atom);
        }
      }
    }
  }

  return features;
}

} // namespace lahuta
