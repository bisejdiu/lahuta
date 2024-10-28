#include "residues.hpp"
#include "GraphMol/MonomerInfo.h"
#include "find_rings.hpp"

namespace lahuta {

std::vector<Residue> Residues::filter(std::function<bool(const Residue &)> predicate) const {
  std::vector<Residue> result;
  for (const auto &residue : residues_) {
    if (predicate(residue)) {
      result.push_back(residue);
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

void Residues::build_residues(const RDKit::RWMol &mol) {
  std::map<std::tuple<std::string, int, std::string>, Residue> residue_map;

  for (const auto &atom : mol.atoms()) {
    if (atom->getAtomicNum() == 1) continue;
    auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    if (!info) continue;

    std::string chain_id = info->getChainId();
    int res_num = info->getResidueNumber();
    std::string res_name = info->getResidueName();

    auto key = std::make_tuple(chain_id, res_num, res_name);

    Residue &residue = residue_map[key];
    if (residue.atoms.empty()) {
      residue = Residue(chain_id, res_num, res_name);
    }
    residue.atoms.push_back(atom);
  }

  residues_.reserve(residue_map.size());
  for (auto &kv : residue_map) {
    residues_.push_back(std::move(kv.second));
  }

  for (const auto &residue : residues_) {
    residues_by_name_[residue.name].push_back(&residue);
  }
}

namespace residue_props {

template <typename ResultType>
std::vector<ResultType> get_aromatic_rings(const Residues &residues, RingProcFunc<ResultType> func) {
  std::vector<ResultType> ring_list;

  // clang-format off
    const std::unordered_map<std::string, std::vector<int>> residue_ring_sizes = {
        {"PHE", {6}},
        {"TYR", {6}},
        {"HIS", {5}},
        {"TRP", {5, 6}}
    };
  // clang-format on

  for (const auto &item : residue_ring_sizes) {
    const std::string &res_name = item.first;
    const std::vector<int> &ring_sizes = item.second;

    auto residues_it = residues.residue_map().find(res_name);
    if (residues_it != residues.residue_map().end()) {
      const std::vector<const Residue *> &res_list = residues_it->second;
      for (const Residue *residue : res_list) {
        for (int ring_size : ring_sizes) {
          ResultType processed_data = func(*residue, ring_size);
          ring_list.push_back(std::move(processed_data));
        }
      }
    }
  }

  return ring_list;
}

template std::vector<std::vector<int>>
get_aromatic_rings<std::vector<int>>(const Residues &, std::function<std::vector<int>(const Residue &, int)>);

template std::vector<std::vector<const RDKit::Atom *>> //
get_aromatic_rings<std::vector<const RDKit::Atom *>>(
    const Residues &, std::function<std::vector<const RDKit::Atom *>(const Residue &, int)>);

/*std::vector<Residue>*/
/*get_unknown_residues(const Residues &residues, const std::set<std::string> &KnownResiduesSet) {*/
/*  return residues.filter([&KnownResiduesSet](const Residue &residue) {*/
/*    return KnownResiduesSet.find(residue.name) == KnownResiduesSet.end();*/
/*  });*/
/*}*/
/**/
/*std::vector<int> // collect all atom indices from residues into one vector*/
/*get_unknown_residues(const Residues &residues, const std::set<std::string> &KnownResiduesSet) {*/
/*  std::vector<int> unknown_indices;*/
/*  for (const auto &residue : get_unknown_residues(residues, KnownResiduesSet)) {*/
/*    for (const auto &atom : residue.atoms) {*/
/*      unknown_indices.push_back(atom->getIdx());*/
/*    }*/
/*  }*/
/*  return unknown_indices;*/
/*}*/

template <typename ReturnType>
ReturnType get_unknown_residues(const Residues &residues, const std::set<std::string> &KnownResiduesSet);

template <>
std::vector<Residue> get_unknown_residues<std::vector<Residue>>(
    const Residues &residues, const std::set<std::string> &KnownResiduesSet) {
  std::vector<Residue> unknown_residues;
  std::copy_if(
      residues.begin(),
      residues.end(),
      std::back_inserter(unknown_residues),
      [&KnownResiduesSet](const Residue &residue) {
        return KnownResiduesSet.find(residue.name) == KnownResiduesSet.end();
      });
  return unknown_residues;
}

template <>
std::vector<int> get_unknown_residues<std::vector<int>>(
    const Residues &residues, const std::set<std::string> &KnownResiduesSet) {
  std::vector<int> unknown_indices;
  for (const auto &residue : get_unknown_residues<std::vector<Residue>>(residues, KnownResiduesSet)) {
    for (const auto &atom : residue.atoms) {
      unknown_indices.push_back(atom->getIdx());
    }
  }
  return unknown_indices;
}

std::vector<std::vector<int>> find_aromatic_rings(const RDKit::RWMol &mol, const Residues &residues) {
  return get_aromatic_rings<std::vector<int>>(residues, [&mol](const Residue &residue, int ring_size) {
    std::vector<int> atom_ids;
    for (const auto &atom : FastRingFinder::find_ring_in_residue(mol, residue, ring_size)) {
      atom_ids.push_back(atom->getIdx());
    }
    return atom_ids;
  });
}

} // namespace residue_props

} // namespace lahuta
