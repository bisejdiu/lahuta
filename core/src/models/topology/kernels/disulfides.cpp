/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string first, last, domain;
 *   std::tie(first, last, domain) = std::make_tuple("besian", "sejdiu", "gmail.com");
 *   return first + last + "@" + domain;
 * }();
 *
 */

#include <rdkit/GraphMol/RWMol.h>

#include "compute/context.hpp"
#include "compute/result.hpp"
#include "models/ssbonds.hpp"
#include "models/topology/data.hpp"
#include "models/topology/kernels.hpp"

// clang-format off
namespace lahuta::models::topology {
namespace C = lahuta::compute;

template <typename DataT>
C::ComputationResult
ModelDisulfidesKernel::execute(C::DataContext<DataT, C::Mut::ReadWrite> &context, const ModelDisulfidesParams &params) {
  auto &data = context.data();

  try {
    if (data.sulphur_atom_indices.empty()) { data.disulfide_pairs.clear(); return C::ComputationResult(true); }

    auto disulfide_pairs = find_disulfide_bonds(data.sulphur_atom_indices, data.conf->getPositions());

    for (const auto &pair : disulfide_pairs) {
      RDKit::Bond *bond = data.bond_pool->createBond(pair.first, pair.second, RDKit::Bond::SINGLE);
      data.bonds.push_back(bond);
    }

    data.disulfide_pairs = std::move(disulfide_pairs);

    return C::ComputationResult(true);
  } catch (const std::exception &e) {
    return C::ComputationResult(C::ComputationError(std::string("Error processing model disulfides: ") + e.what()));
  }
}

template C::ComputationResult ModelDisulfidesKernel::execute<ModelData>(C::DataContext<ModelData, C::Mut::ReadWrite> &, const ModelDisulfidesParams &);

} // namespace lahuta::models::topology
