#include "compute/context.hpp"
#include "compute/result.hpp"
#include "models/ssbonds.hpp"
#include "models/topology/data.hpp"
#include "models/topology/kernels.hpp"
#include <rdkit/GraphMol/RWMol.h>

// clang-format off
namespace lahuta::models::topology {

template <typename DataT>
ComputationResult
ModelDisulfidesKernel::execute(DataContext<DataT, Mut::ReadWrite> &context, const ModelDisulfidesParams &params) {
  auto &data = context.data();

  try {
    if (data.sulphur_atom_indices.empty()) { data.disulfide_pairs.clear(); return ComputationResult(true); }

    auto disulfide_pairs = find_disulfide_bonds(data.sulphur_atom_indices, data.conf->getPositions());

    for (const auto &pair : disulfide_pairs) {
      RDKit::Bond *bond = data.bond_pool->createBond(pair.first, pair.second, RDKit::Bond::SINGLE);
      data.bonds.push_back(bond);
    }

    data.disulfide_pairs = std::move(disulfide_pairs);

    return ComputationResult(true);
  } catch (const std::exception &e) {
    return ComputationResult(ComputationError(std::string("Error processing model disulfides: ") + e.what()));
  }
}

template ComputationResult ModelDisulfidesKernel::execute<ModelData>(DataContext<ModelData, Mut::ReadWrite> &, const ModelDisulfidesParams &);

} // namespace lahuta::models::topology
