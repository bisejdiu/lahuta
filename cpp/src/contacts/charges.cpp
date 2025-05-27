#include "contacts/charges.hpp"
#include "contacts/heuristics.hpp"
#include "contacts/utils.hpp"
#include "entities/records.hpp"
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

std::vector<GroupRec> add_positive_charges(const RDKit::RWMol &mol, const Residues &residues) {
  std::optional<std::unordered_map<const RDKit::Atom *, FeatureGroup>> groups;
  std::unordered_set<const RDKit::Atom *> added_atoms;

  std::vector<GroupRec> features;
  for (const auto &residue : residues) {
    // Handle positively charged residues (ARG, HIS, LYS)
    if (definitions::is_positive_charge(residue.name)) {
      std::vector<uint32_t> member_indices;
      for (const auto *atom : residue.atoms) {
        auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
        std::string atom_name = res_info->getName();

        if (atom->getAtomicNum() == Element::N && definitions::ProteinBackboneAtoms.count(atom_name) == 0) {
          member_indices.push_back(static_cast<uint32_t>(atom->getIdx()));
          added_atoms.insert(atom);
        }
      }

      if (!member_indices.empty()) {
        features.push_back(GroupRec{
          /*.a_type =*/ AtomType::POS_IONISABLE,
          /*.type   =*/ FeatureGroup::None,
          /*.atoms  =*/ std::move(member_indices),
          /*.center =*/ RDGeom::Point3D(0,0,0)
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

          std::vector<uint32_t> nitrogen_indices;
          nitrogen_indices.reserve(nitrogens.size());
          for (const auto* n_atom : nitrogens) {
            nitrogen_indices.push_back(static_cast<uint32_t>(n_atom->getIdx()));
            added_atoms.insert(n_atom);
          }

          features.push_back(GroupRec{
            /*.a_type =*/ AtomType::POS_IONISABLE,
            /*.type   =*/ it->second,
            /*.atoms  =*/ std::move(nitrogen_indices),
            /*.center =*/ RDGeom::Point3D(0,0,0)
          });
        } else {
          // Add remaining positively charged atoms not already added
          if (atom->getFormalCharge() > 0 && added_atoms.count(atom) == 0) {
            std::vector<uint32_t> atom_indices = {static_cast<uint32_t>(atom->getIdx())};
            features.push_back(GroupRec{
              /*.a_type =*/ AtomType::POS_IONISABLE,
              /*.type   =*/ FeatureGroup::None,
              /*.atoms  =*/ std::move(atom_indices),
              /*.center =*/ RDGeom::Point3D(0,0,0)
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
      std::vector<uint32_t> member_indices;
      for (const auto *atom : residue.atoms) {
        auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
        std::string atom_name = res_info->getName();

        if (atom->getAtomicNum() == Element::O && definitions::ProteinBackboneAtoms.count(atom_name) == 0) {
          member_indices.push_back(static_cast<uint32_t>(atom->getIdx()));
          added_atoms.insert(atom);
        }
      }

      if (member_indices.empty()) continue;
      
      features.push_back(GroupRec{
        /*.a_type =*/ AtomType::NEG_IONISABLE,
        /*.type   =*/ FeatureGroup::None,
        /*.atoms  =*/ std::move(member_indices),
        /*.center =*/ RDGeom::Point3D(0,0,0)
      });
    } else if (definitions::is_base(residue.name)) {
      // Handle nucleic acid bases (DNA/RNA)
      for (const auto *atom : residue.atoms) {
        if (is_P_in_phosphate(mol, *atom)) {
          auto oxygens = bonded_atoms(mol, atom, Element::O);
          if (oxygens.empty()) continue;

          std::vector<uint32_t> oxygen_indices;
          oxygen_indices.reserve(oxygens.size());
          for (const auto* o_atom : oxygens) {
            oxygen_indices.push_back(static_cast<uint32_t>(o_atom->getIdx()));
            added_atoms.insert(o_atom);
          }

          features.push_back(GroupRec{
            /*.a_type =*/ AtomType::NEG_IONISABLE,
            /*.type   =*/ FeatureGroup::Phosphate,
            /*.atoms  =*/ std::move(oxygen_indices),
            /*.center =*/ RDGeom::Point3D(0,0,0)
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

          std::vector<uint32_t> oxygen_indices;
          oxygen_indices.reserve(oxygens.size());
          for (const auto* o_atom : oxygens) {
            oxygen_indices.push_back(static_cast<uint32_t>(o_atom->getIdx()));
            added_atoms.insert(o_atom);
          }

          features.push_back(GroupRec{
            /*.a_type =*/ AtomType::NEG_IONISABLE,
            /*.type   =*/ it->second,
            /*.atoms  =*/ std::move(oxygen_indices),
            /*.center =*/ RDGeom::Point3D(0,0,0)
          });
        } else {
          // Collect non-backbone nitrogen atoms
          auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
          std::string atom_name = res_info->getName();

          if (atom->getAtomicNum() == Element::N && definitions::ProteinBackboneAtoms.count(atom_name) == 0) {
            std::vector<uint32_t> atom_indices = {static_cast<uint32_t>(atom->getIdx())};
            features.push_back(GroupRec{
              /*.a_type =*/ AtomType::NEG_IONISABLE,
              /*.type   =*/ FeatureGroup::None,
              /*.atoms  =*/ std::move(atom_indices),
              /*.center =*/ RDGeom::Point3D(0,0,0)
            });
            added_atoms.insert(atom);
          }
        }

        // Add remaining negatively charged atoms not already added
        if (atom->getFormalCharge() < 0 && added_atoms.count(atom) == 0) {
          std::vector<uint32_t> atom_indices = {static_cast<uint32_t>(atom->getIdx())};
          features.push_back(GroupRec{
            /*.a_type =*/ AtomType::NEG_IONISABLE,
            /*.type   =*/ FeatureGroup::None,
            /*.atoms  =*/ std::move(atom_indices),
            /*.center =*/ RDGeom::Point3D(0,0,0)
          });
          added_atoms.insert(atom);
        }
      }
    }
  }

  return features;
}

} // namespace lahuta
