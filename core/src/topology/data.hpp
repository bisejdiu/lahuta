#pragma once

#include "entities/records.hpp"
#include "residues.hpp"
#include <memory>
#include <vector>
#include <rdkit/GraphMol/RWMol.h>

// clang-format off
namespace lahuta::topology {

namespace {
template<typename T>
struct always_false : std::false_type {};
}

struct TopologyData {
  TopologyData(std::shared_ptr<RDKit::RWMol> mol_ptr)
      : mol(mol_ptr), residues(std::make_unique<Residues>(*mol)) {}

  template <typename T>
  const std::vector<T> &get_container() const {

    if constexpr (std::is_same_v<T, AtomRec>) { return atoms; }
    else if constexpr (std::is_same_v<T, RingRec>) { return rings; }
    else if constexpr (std::is_same_v<T, GroupRec>) { return groups; }
    else {
      static_assert(always_false<T>::value, "Unsupported record type");
    }
  }

  std::shared_ptr<RDKit::RWMol> mol;
  std::unique_ptr<Residues> residues;
  std::vector<AtomRec>  atoms;
  std::vector<RingRec>  rings;
  std::vector<GroupRec> groups;
};

} // namespace lahuta::topology
