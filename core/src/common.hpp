#ifndef LAHUTA_COMMON_HPP
#define LAHUTA_COMMON_HPP

#include "GraphMol/MonomerInfo.h"
#include "GraphMol/RWMol.h"
#include "logging.hpp"

namespace lahuta {

struct AtomInfo {
  RDKit::Atom *atom;
  const RDKit::AtomPDBResidueInfo *info;
  bool is_hydrogen;

  AtomInfo(RDKit::Atom *atom, const RDKit::AtomPDBResidueInfo *info, bool is_hydrogen)
      : atom(atom), info(info), is_hydrogen(is_hydrogen) {}

  AtomInfo(const RDKit::RWMol &mol, int idx)
      : atom(const_cast<RDKit::Atom *>(mol.getAtomWithIdx(idx))),
        info(static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo())),
        is_hydrogen(atom->getAtomicNum() == 1) {}
};

inline auto is_same_conformer(const AtomInfo &a, const AtomInfo &b) {
  return a.info->getAltLoc().empty() || b.info->getAltLoc().empty()
         || a.info->getAltLoc() == b.info->getAltLoc();
};

namespace common {

/// Check if a container contains a value
template <typename Container, typename T> bool contains(const Container &container, const T &value) {
  return std::find(std::begin(container), std::end(container), value) != std::end(container);
}

struct PairHash {
  std::size_t operator()(const std::pair<int, int> &p) const {
    return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
  }
};

inline bool is_ring_aromatic(const RDKit::RWMol &mol, const RDKit::INT_VECT &ring) {
  return std::all_of(ring.begin(), ring.end(), [&mol](int idx) {
    return mol.getAtomWithIdx(idx)->getIsAromatic();
  });
}

inline bool has_any_aromatic_atom(const RDKit::RWMol &mol, const RDKit::INT_VECT &ring) {
  return std::any_of(ring.begin(), ring.end(), [&mol](int idx) {
    return mol.getAtomWithIdx(idx)->getIsAromatic();
  });
}

inline std::vector<int> factorize(const std::vector<std::string> &labels) {
  std::vector<int> ids(labels.size());

  // hash map from labels to ids
  std::unordered_map<std::string_view, int> label_to_id;
  label_to_id.reserve(labels.size());

  int current_id = 0;
  for (size_t i = 0; i < labels.size(); ++i) {
    std::string_view label = labels[i];
    auto it = label_to_id.find(label);
    if (it == label_to_id.end()) {
      label_to_id[label] = current_id;
      ids[i] = current_id;
      ++current_id;
    } else {
      ids[i] = it->second;
    }
  }

  return ids;
}

inline int count_unique(const std::vector<int> &vec) {
  std::unordered_set<int> unique_elements(vec.begin(), vec.end());
  return unique_elements.size();
}

inline int count_unique(const std::vector<std::string> &vec) {
  std::unordered_set<std::string_view> unique_elements;
  unique_elements.reserve(vec.size());

  for (const auto &str : vec) {
    unique_elements.insert(std::string_view(str));
  }

  return unique_elements.size();
}

inline void log_atom_info(const RDKit::Atom *atom) {
  auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
  if (info) {
    Logger::get_logger()->info(
        "Atom: {} {} {} {} {}",
        atom->getIdx(),
        info->getName(),
        info->getResidueName(),
        info->getResidueNumber(),
        info->getChainId());
  } else {
    Logger::get_logger()->trace("Atom: {} {}", atom->getIdx(), atom->getSymbol());
  }
}

inline void log_bond_info(const RDKit::Bond *bond) {
  auto *begin = bond->getBeginAtom();
  auto *end = bond->getEndAtom();
  auto *info_begin = static_cast<const RDKit::AtomPDBResidueInfo *>(begin->getMonomerInfo());
  auto *info_end = static_cast<const RDKit::AtomPDBResidueInfo *>(end->getMonomerInfo());
  if (info_begin && info_end) {
    Logger::get_logger()->info(
        "Bond: {}-{}-{}-{}-{} {}-{}-{}-{}-{}",
        begin->getIdx(),
        info_begin->getName(),
        info_begin->getResidueName(),
        info_begin->getResidueNumber(),
        info_begin->getChainId(),
        end->getIdx(),
        info_end->getName(),
        info_end->getResidueName(),
        info_end->getResidueNumber(),
        info_end->getChainId());
  }
}

inline void log_bond_info(const RDKit::Atom *a1, const RDKit::Atom *a2) {
  auto *info_a1 = static_cast<const RDKit::AtomPDBResidueInfo *>(a1->getMonomerInfo());
  auto *info_a2 = static_cast<const RDKit::AtomPDBResidueInfo *>(a2->getMonomerInfo());
  if (info_a1 && info_a2) {
    Logger::get_logger()->info(
        "Bond: {}-{}-{}-{}-{} {}-{}-{}-{}-{}",
        a1->getIdx(),
        info_a1->getName(),
        info_a1->getResidueName(),
        info_a1->getResidueNumber(),
        info_a1->getChainId(),
        a2->getIdx(),
        info_a2->getName(),
        info_a2->getResidueName(),
        info_a2->getResidueNumber(),
        info_a2->getChainId());
  }
}

} // namespace common
} // namespace lahuta

#endif // LAHUTA_COMMON_HPP
