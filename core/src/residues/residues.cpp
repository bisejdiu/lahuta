/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: std::string{"besian"} + "sejdiu" + "@gmail.com";
 *
 */

#include <rdkit/GraphMol/MonomerInfo.h>

#include "hash/fnv1a.hpp"
#include "residues/definitions.hpp"
#include "residues/residues.hpp"
#include "rings/find_rings.hpp"

namespace lahuta {

namespace {

struct ResidueKey {
  std::uint64_t key;

  ResidueKey(std::string_view chain_id, int res_num, std::string_view res_name) noexcept {
    const std::uint32_t chain_hash = fnv1a_32(chain_id);
    const std::uint32_t name_hash  = fnv1a_32(res_name);

    key = (static_cast<std::uint64_t>(chain_hash & 0xFFFF) << 48) |
          (static_cast<std::uint64_t>(static_cast<std::uint32_t>(res_num)) << 16) |
          (static_cast<std::uint64_t>(name_hash & 0xFFFF));
  }

  bool operator==(const ResidueKey &other) const noexcept { return key == other.key; }
};

struct ResidueKeyHash {
  std::size_t operator()(const ResidueKey &k) const noexcept { return static_cast<std::size_t>(k.key); }
};

// Originals are stored to handle fallback in case of collisions
struct ResidueEntry {
  std::size_t index;
  std::string chain_id;
  int res_num;
  std::string res_name;
};

} // namespace

bool Residues::build() {
  bool success = false;
  build_residues(mol_, success);
  return success;
}

Residues Residues::filter(std::function<bool(const Residue &)> predicate) const {
  Residues result(mol_);
  result.atom_to_residue_idx_.assign(atom_to_residue_idx_.size(), -1);
  for (const auto &residue : residues_) {
    if (predicate(residue)) {
      const int mapped_idx = static_cast<int>(result.residues_.size());
      result.residues_.push_back(residue);
      for (const auto *atom : residue.atoms) {
        result.atom_to_residue_idx_[static_cast<std::size_t>(atom->getIdx())] = mapped_idx;
      }
    }
  }
  return result;
}

Residues Residues::filter(std::function<bool(const std::string &)> predicate) const {
  Residues result(mol_);
  result.atom_to_residue_idx_.assign(atom_to_residue_idx_.size(), -1);
  for (const auto &residue : residues_) {
    if (predicate(residue.name)) {
      const int mapped_idx = static_cast<int>(result.residues_.size());
      result.residues_.push_back(residue);
      for (const auto *atom : residue.atoms) {
        result.atom_to_residue_idx_[static_cast<std::size_t>(atom->getIdx())] = mapped_idx;
      }
    }
  }
  return result;
}

template <typename ResultType>
std::vector<ResultType> Residues::map(std::function<ResultType(const Residue &)> func) const {
  std::vector<ResultType> result;
  result.reserve(residues_.size());
  for (const auto &residue : residues_) {
    result.push_back(func(residue));
  }
  return result;
}

void Residues::build_residues(const RDKit::RWMol &mol, bool &status) {
  const std::size_t num_atoms = mol.getNumAtoms();
  if (num_atoms == 0) {
    residues_.clear();
    atom_to_residue_idx_.clear();
    status = true;
    return;
  }

  const std::size_t estimated_residues = (num_atoms / 15) + 1; // ~15 atoms per residue

  std::vector<Residue> residues;
  residues.reserve(estimated_residues);
  std::vector<int> atom_to_residue_idx(num_atoms, -1);

  std::unordered_map<ResidueKey, ResidueEntry, ResidueKeyHash> residue_map;
  residue_map.reserve(estimated_residues);

  for (const auto &atom : mol.atoms()) {
    auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    if (!info) continue;

    const int res_num               = info->getResidueNumber();
    const std::string &res_name_ref = info->getResidueName();
    const std::string &chain_id_ref = info->getChainId();

    const std::string_view res_name_sv(res_name_ref);
    const std::string_view chain_id_sv(chain_id_ref);

    const ResidueKey key(chain_id_sv, res_num, res_name_sv);
    const std::size_t atom_idx = static_cast<std::size_t>(atom->getIdx());
    auto it = residue_map.find(key);

    if (it != residue_map.end()) {
      auto &entry = it->second;
      if (entry.chain_id == chain_id_sv && entry.res_num == res_num && entry.res_name == res_name_sv) {
        residues[entry.index].atoms.push_back(atom);
        atom_to_residue_idx[atom_idx] = static_cast<int>(entry.index);
      } else {
        bool found = false;
        for (std::size_t i = 0; i < residues.size(); ++i) {
          auto &r = residues[i];
          if (r.number == res_num && r.chain_id == chain_id_sv && r.name == res_name_sv) {
            r.atoms.push_back(atom);
            atom_to_residue_idx[atom_idx] = static_cast<int>(i);
            found = true;
            break;
          }
        }
        if (!found) {
          const std::size_t new_idx = residues.size();
          residues.emplace_back(std::string(chain_id_sv), res_num, std::string(res_name_sv), "");
          residues.back().atoms.push_back(atom);
          atom_to_residue_idx[atom_idx] = static_cast<int>(new_idx);
        }
      }
    } else {
      const std::size_t idx = residues.size();
      residues.emplace_back(std::string(chain_id_sv), res_num, std::string(res_name_sv), "");
      residues.back().atoms.push_back(atom);
      atom_to_residue_idx[atom_idx] = static_cast<int>(idx);
      residue_map.emplace(key,
                          ResidueEntry{idx, std::string(chain_id_sv), res_num, std::string(res_name_sv)});
    }
  }

  for (std::size_t i = 0; i < residues.size(); ++i) {
    residues[i].idx = static_cast<unsigned>(i);
  }

  residues_ = std::move(residues);
  atom_to_residue_idx_ = std::move(atom_to_residue_idx);
  status    = true;
}

std::vector<std::vector<int>> find_and_process_aromatic_residues(const RDKit::RWMol &mol,
                                                                 const Residues &residues) {

  using RingSize               = definitions::arom_rings::RingSize;
  constexpr auto &AromRingSize = definitions::arom_rings::AromaticResiduesRingSizes;

  std::vector<std::vector<int>> ring_vector;
  for (const auto &residue : residues) {

    // clang-format off
    RingSize ring_size;
    if      (residue.name == "PHE") { ring_size = RingSize::RS_6; }
    else if (residue.name == "TYR") { ring_size = RingSize::RS_6; }
    else if (residue.name == "HIS") { ring_size = RingSize::RS_5; }
    else if (residue.name == "TRP") { ring_size = static_cast<RingSize>(RingSize::RS_5 | RingSize::RS_6); }
    // clang-format on
    else {
      // linear search through remaining predefined residues
      auto it = std::find_if(AromRingSize.begin(), AromRingSize.end(), [&](const auto &e) {
        return e.first == residue.name;
      });
      if (it == std::end(AromRingSize)) continue;

      ring_size = it->second;
    }

    std::vector<int> ring_sizes = definitions::arom_rings::get_ringsizes(ring_size);
    for (int ring_size : ring_sizes) {
      if (residue.atoms.size() < static_cast<size_t>(ring_size)) continue;

      std::vector<int> processed_data = FastRingFinder::find_ring_in_residue(residues.get_mol(),
                                                                             residue,
                                                                             ring_size);
      if (processed_data.size() != ring_size) continue;

      ring_vector.push_back(std::move(processed_data));
    }
  }

  return ring_vector;
}

} // namespace lahuta
