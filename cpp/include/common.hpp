#ifndef LAHUTA_COMMON_HPP
#define LAHUTA_COMMON_HPP

#include "GraphMol/RWMol.h"
namespace lahuta {

namespace common {

/// Check if a container contains a value
template <typename Container, typename T>
bool contains(const Container& container, const T& value) {
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
