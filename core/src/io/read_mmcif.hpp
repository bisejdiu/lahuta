/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   using Part = std::variant<const char*, std::string_view>;
 *   std::array<Part, 3> parts{Part{"besian"}, Part{"sejdiu"}, Part{"@gmail.com"}};
 *   std::string s;
 *   for (const auto& p : parts) {
 *     std::visit([&s](auto&& arg) { s += arg; }, p);
 *   }
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_READ_MMCIF_HPP
#define LAHUTA_READ_MMCIF_HPP

#include <cctype>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <gemmi/atox.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/elem.hpp>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>

#include "gemmi/numb.hpp"
#include "logging/logging.hpp"

namespace lahuta {
namespace detail {

struct ResidueAtomRef {
  unsigned int atom_idx = 0;
  char alt_loc          = '\0';
};

struct ResidueAtomBucket {
  std::unordered_map<std::string, std::vector<ResidueAtomRef>> atoms_by_name;
};

inline RDKit::AtomMonomerInfo &require_monomer_info(RDKit::Atom &atom) {
  auto *info = atom.getMonomerInfo();
  if (!info) {
    throw std::runtime_error("mmCIF component chemistry requires monomer info on every atom");
  }
  return *info;
}

inline void trim_mmcif_field(std::string &value) noexcept {
  if (value.empty()) return;

  const unsigned char first = static_cast<unsigned char>(value.front());
  const unsigned char last  = static_cast<unsigned char>(value.back());
  if (first > ' ' && last > ' ') {
    return;
  }

  std::size_t start = 0;
  std::size_t end   = value.size();

  while (start < end && static_cast<unsigned char>(value[start]) <= ' ') {
    ++start;
  }
  while (end > start && static_cast<unsigned char>(value[end - 1]) <= ' ') {
    --end;
  }

  if (start == 0 && end == value.size()) {
    return;
  }
  if (start == end) {
    value.clear();
    return;
  }
  if (end < value.size()) {
    value.erase(end);
  }
  if (start > 0) {
    value.erase(0, start);
  }
}

inline bool mmcif_yes(std::string value) noexcept {
  trim_mmcif_field(value);
  if (value.empty()) return false;
  const char c = static_cast<char>(std::toupper(static_cast<unsigned char>(value.front())));
  return c == 'Y';
}

inline RDKit::Bond::BondType mmcif_bond_type(std::string value) noexcept {
  trim_mmcif_field(value);
  if (value == "doub") return RDKit::Bond::BondType::DOUBLE;
  if (value == "trip") return RDKit::Bond::BondType::TRIPLE;
  if (value == "quad") return RDKit::Bond::BondType::QUADRUPLE;
  if (value == "arom") return RDKit::Bond::BondType::AROMATIC;
  if (value == "delo") return RDKit::Bond::BondType::AROMATIC;
  return RDKit::Bond::BondType::SINGLE;
}

inline bool mmcif_alt_compatible(char a, char b) noexcept { return a == '\0' || b == '\0' || a == b; }

inline std::string make_residue_instance_key(std::string_view comp_id, std::string_view chain_id,
                                             std::string_view auth_seq_id, std::string_view ins_code) {
  std::string key;
  key.reserve(comp_id.size() + chain_id.size() + auth_seq_id.size() + ins_code.size() + 4);
  key.append(comp_id);
  key.push_back('\x1f');
  key.append(chain_id);
  key.push_back('\x1f');
  key.append(auth_seq_id);
  key.push_back('\x1f');
  key.append(ins_code);
  return key;
}

inline void
apply_mmcif_component_aromaticity(RDKit::RWMol &mol, gemmi::cif::Block &block,
                                  const std::unordered_map<std::string, ResidueAtomBucket> &residue_atoms) {
  namespace cif = gemmi::cif;

  cif::Table atom_table = block.find("_chem_comp_atom.", {"comp_id", "atom_id", "?pdbx_aromatic_flag"});
  if (atom_table.length() == 0) return;

  enum { kCompId = 0, kAtomId, kAromaticFlag };

  for (auto row : atom_table) {
    if (!row.has2(kCompId) || !row.has2(kAtomId)) continue;

    std::string comp_id = row.str(kCompId);
    std::string atom_id = row.str(kAtomId);
    detail::trim_mmcif_field(comp_id);
    detail::trim_mmcif_field(atom_id);
    if (comp_id.empty() || atom_id.empty()) continue;
    if (!row.has2(kAromaticFlag) || !detail::mmcif_yes(row.str(kAromaticFlag))) continue;

    const std::string prefix = comp_id + '\x1f';
    for (const auto &[instance_key, bucket] : residue_atoms) {
      if (instance_key.rfind(prefix, 0) != 0) continue;
      auto it = bucket.atoms_by_name.find(atom_id);
      if (it == bucket.atoms_by_name.end()) continue;
      for (const auto &ref : it->second) {
        mol.getAtomWithIdx(ref.atom_idx)->setIsAromatic(true);
      }
    }
  }
}

inline void
apply_mmcif_component_bonds(RDKit::RWMol &mol, gemmi::cif::Block &block,
                            const std::unordered_map<std::string, ResidueAtomBucket> &residue_atoms) {
  namespace cif = gemmi::cif;

  cif::Table bond_table = block.find(
      "_chem_comp_bond.",
      {"comp_id", "atom_id_1", "atom_id_2", "?value_order", "?pdbx_aromatic_flag"});
  if (bond_table.length() == 0) return;

  enum { kCompId = 0, kAtomId1, kAtomId2, kValueOrder, kAromaticFlag };

  for (auto row : bond_table) {
    if (!row.has2(kCompId) || !row.has2(kAtomId1) || !row.has2(kAtomId2)) continue;

    std::string comp_id   = row.str(kCompId);
    std::string atom_id_1 = row.str(kAtomId1);
    std::string atom_id_2 = row.str(kAtomId2);
    detail::trim_mmcif_field(comp_id);
    detail::trim_mmcif_field(atom_id_1);
    detail::trim_mmcif_field(atom_id_2);
    if (comp_id.empty() || atom_id_1.empty() || atom_id_2.empty()) continue;

    const auto bond_type = row.has2(kValueOrder) ? detail::mmcif_bond_type(row.str(kValueOrder))
                                                 : RDKit::Bond::BondType::SINGLE;
    const bool aromatic  = bond_type == RDKit::Bond::BondType::AROMATIC ||
                          (row.has2(kAromaticFlag) && detail::mmcif_yes(row.str(kAromaticFlag)));

    const std::string prefix = comp_id + '\x1f';
    for (const auto &[instance_key, bucket] : residue_atoms) {
      if (instance_key.rfind(prefix, 0) != 0) continue;

      for (const auto &[_, refs] : bucket.atoms_by_name) {
        for (const auto &ref : refs) {
          auto *atom = mol.getAtomWithIdx(ref.atom_idx);
          detail::require_monomer_info(*atom).setHasChemCompBondSchema(true);
        }
      }

      auto it_a = bucket.atoms_by_name.find(atom_id_1);
      auto it_b = bucket.atoms_by_name.find(atom_id_2);
      if (it_a == bucket.atoms_by_name.end() || it_b == bucket.atoms_by_name.end()) continue;

      for (const auto &a_ref : it_a->second) {
        for (const auto &b_ref : it_b->second) {
          if (a_ref.atom_idx == b_ref.atom_idx) continue;
          if (!detail::mmcif_alt_compatible(a_ref.alt_loc, b_ref.alt_loc)) continue;

          RDKit::Bond *bond = mol.getBondBetweenAtoms(a_ref.atom_idx, b_ref.atom_idx);
          if (!bond) {
            mol.addBond(a_ref.atom_idx, b_ref.atom_idx, bond_type);
            bond = mol.getBondBetweenAtoms(a_ref.atom_idx, b_ref.atom_idx);
          } else {
            bond->setBondType(bond_type);
          }

          if (bond) {
            bond->setIsAromatic(aromatic);
          }
        }
      }
    }
  }
}

} // namespace detail

inline std::shared_ptr<RDKit::RWMol> make_mol_from_block(const gemmi::cif::Block &block_,
                                                         const std::string &source) {
  namespace cif = gemmi::cif;

  auto &block = const_cast<cif::Block &>(block_);
  auto mol    = std::make_shared<RDKit::RWMol>();

  enum {
    kId = 0,
    kGroupPdb,
    kSymbol,
    kLabelAtomId,
    kAltId,
    kLabelCompId,
    kLabelAsymId,
    kLabelEntityId,
    kLabelSeqId,
    kInsCode,
    kX,
    kY,
    kZ,
    kOcc,
    kBiso,
    kCharge,
    kAuthSeqId,
    kAuthCompId,
    kAuthAsymId,
    kAuthAtomId,
    kModelNum
  };

  cif::Table atom_table = block.find("_atom_site.",
                                     {"id",
                                      "?group_PDB",
                                      "type_symbol",
                                      "?label_atom_id",
                                      "label_alt_id",
                                      "?label_comp_id",
                                      "label_asym_id",
                                      "?label_entity_id",
                                      "?label_seq_id",
                                      "?pdbx_PDB_ins_code",
                                      "Cartn_x",
                                      "Cartn_y",
                                      "Cartn_z",
                                      "occupancy",
                                      "B_iso_or_equiv",
                                      "?pdbx_formal_charge",
                                      "auth_seq_id",
                                      "?auth_comp_id",
                                      "?auth_asym_id",
                                      "?auth_atom_id",
                                      "?pdbx_PDB_model_num"});

  if (atom_table.length() == 0) throw std::runtime_error("No _atom_site table found in the mmCIF file");

  const int kAsymId = atom_table.first_of(kAuthAsymId, kLabelAsymId);
  const int kCompId = atom_table.first_of(kAuthCompId, kLabelCompId);
  const int kAtomId = atom_table.first_of(kAuthAtomId, kLabelAtomId);

  if (!atom_table.has_column(kCompId))
    throw std::runtime_error("Neither _atom_site.label_comp_id nor auth_comp_id found");
  if (!atom_table.has_column(kAtomId))
    throw std::runtime_error("Neither _atom_site.label_atom_id nor auth_atom_id found");

  // we expect that models (if multiple) are ordered by model number (1, 1, 2, 2, 3, 3)
  std::vector<std::vector<RDGeom::Point3D>> model_coords;
  model_coords.reserve(20);
  std::unordered_map<std::string, detail::ResidueAtomBucket> residue_atoms;

  std::string current_model;
  for (auto row : atom_table) {
    std::string model_id = row.has(kModelNum) ? row.str(kModelNum) : "1";

    if (model_coords.empty() || model_id != current_model) {
      current_model = model_id;
      model_coords.emplace_back();
      model_coords.back().reserve(atom_table.length());
    }

    model_coords.back().emplace_back(cif::as_number(row[kX]),
                                     cif::as_number(row[kY]),
                                     cif::as_number(row[kZ]));

    if (model_coords.size() > 1) continue;

    std::string symbol = cif::as_string(row[kSymbol]);
    detail::trim_mmcif_field(symbol);
    int atomic_number = gemmi::Element(symbol).atomic_number();
    auto *rdkit_atom  = new RDKit::Atom(atomic_number);

    if (row.has2(kCharge)) rdkit_atom->setFormalCharge(cif::as_int(row[kCharge]));

    mol->addAtom(rdkit_atom, true, true);

    auto atom_name  = cif::as_string(row[kAtomId]);
    auto res_name   = cif::as_string(row[kCompId]);
    auto chain_name = cif::as_string(row[kAsymId]);
    detail::trim_mmcif_field(atom_name);
    detail::trim_mmcif_field(res_name);
    detail::trim_mmcif_field(chain_name);
    auto auth_seq_id = cif::as_string(row[kAuthSeqId]);
    int seq_id       = gemmi::string_to_int(auth_seq_id, false);
    char alt_loc     = cif::as_char(row[kAltId], '\0');
    auto alt_loc_str = (alt_loc == '\0') ? "" : std::string(1, alt_loc);
    int serial       = gemmi::string_to_int(row[kId], false);

    auto *info = new RDKit::AtomPDBResidueInfo(atom_name, serial, alt_loc_str, res_name, seq_id, chain_name);

    if (row.has2(kGroupPdb)) {
      char het_flag = '\0';
      for (int i = 0; i < 2; ++i) {
        const char c = toupper(row[kGroupPdb][i]);
        if (c == 'A' || c == 'H') het_flag = c;
      }
      info->setIsHeteroAtom(het_flag == 'H'); // it is false by default
    }

    if (row.has2(kLabelSeqId)) {
      int label_seq = cif::as_int(row[kLabelSeqId]);
      info->setResidueNumber(label_seq);
    }

    info->setMonomerType(RDKit::AtomMonomerInfo::PDBRESIDUE);
    rdkit_atom->setMonomerInfo(info);

    const auto residue_key = detail::make_residue_instance_key(
        res_name,
        chain_name,
        auth_seq_id,
        row.has2(kInsCode) ? std::string(row.str(kInsCode)) : std::string{});
    residue_atoms[residue_key].atoms_by_name[atom_name].push_back(
        detail::ResidueAtomRef{rdkit_atom->getIdx(), alt_loc});
  }

  detail::apply_mmcif_component_aromaticity(*mol, block, residue_atoms);
  detail::apply_mmcif_component_bonds(*mol, block, residue_atoms);

  const std::size_t atoms_per_model = model_coords[0].size();

  auto new_end = std::remove_if(
      model_coords.begin(),
      model_coords.end(),
      [atoms_per_model](auto const &coords) { return coords.size() != atoms_per_model; });

  if (new_end != model_coords.end()) {
    model_coords.erase(new_end, model_coords.end());
    Logger::get_logger()->warn("Dropping {} models with inconsistent number of atoms from '{}'",
                               model_coords.size(),
                               source);
  }

  unsigned int conf_id = 0;
  for (auto &coords : model_coords) {
    RDKit::Conformer *conformer = new RDKit::Conformer(atoms_per_model);
    conformer->setId(conf_id++);
    conformer->setAllAtomPositions(std::move(coords));
    mol->addConformer(conformer, true);
  }

  if (mol->getNumConformers() != 0) {
    Logger::get_logger()->debug("Loaded {} models with {} atoms",
                                mol->getNumConformers(),
                                mol->getNumAtoms());
  }

  mol->updatePropertyCache(false);
  return mol;
}

} // namespace lahuta

#endif // LAHUTA_READ_MMCIF_HPP
