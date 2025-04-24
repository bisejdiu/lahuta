#pragma once

#include "entities.hpp"
#include "residues.hpp"
#include <memory>
#include <rdkit/GraphMol/RWMol.h>

// clang-format off
namespace lahuta::topology {

struct TopologyData {

  TopologyData(std::shared_ptr<RDKit::RWMol> mol_ptr)
      : mol(mol_ptr), residues(std::make_unique<Residues>(*mol)) {}

  std::shared_ptr<RDKit::RWMol> mol;
  std::unique_ptr<Residues> residues;

  AtomEntityCollection  atom_types;
  RingEntityCollection  rings;
  GroupEntityCollection features;
};

} // namespace lahuta::topology
