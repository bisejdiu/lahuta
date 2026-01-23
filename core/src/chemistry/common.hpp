#ifndef LAHUTA_COMMON_HPP
#define LAHUTA_COMMON_HPP

#include <algorithm>
#include <iterator>

#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>

#include "logging/logging.hpp"

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
