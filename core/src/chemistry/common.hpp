#ifndef LAHUTA_COMMON_HPP
#define LAHUTA_COMMON_HPP

#include <algorithm>
#include <iterator>

#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>

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

} // namespace common
} // namespace lahuta

#endif // LAHUTA_COMMON_HPP
