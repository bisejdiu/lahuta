#ifndef LAHUTA_MODELS_TOPOLOGY_DATA_HPP
#define LAHUTA_MODELS_TOPOLOGY_DATA_HPP

#include "models/factory.hpp"
#include "models/parser.hpp"
#include "models/pools.hpp"
#include <memory>
#include <rdkit/GraphMol/GraphDefs.hpp>
#include <rdkit/GraphMol/RWMol.h>
#include <vector>

// clang-format off
namespace lahuta::models::topology {

struct ModelData {
  ModelData(const ModelParserResult& input)
    : input_data(&input),
      info_pool(PoolFactory<InfoPool>::getFreshPoolForCurrentThread()),
      atom_pool(PoolFactory<AtomPool>::getFreshPoolForCurrentThread()),
      bond_pool(PoolFactory<BondPool>::getFreshPoolForCurrentThread()) {}

  const ModelParserResult* input_data;

  InfoPool* info_pool;
  AtomPool* atom_pool;
  BondPool* bond_pool;

  std::shared_ptr<RDKit::RWMol> mol;
  std::unique_ptr<RDKit::Conformer> conf;

  std::vector<RDKit::Atom*> vertices;
  std::vector<RDKit::Bond*> bonds;
  std::vector<int> sulphur_atom_indices;
  std::vector<std::vector<int>> aromatic_atom_indices;
  std::vector<std::vector<int>> aromatic_bond_indices;
  std::vector<std::pair<int,int>> disulfide_pairs;

  int atom_idx = 0;
  RDKit::GraphType graph_type = RDKit::GraphType::MolGraph;
};

} // namespace lahuta::models::topology

#endif // LAHUTA_MODELS_TOPOLOGY_DATA_HPP
