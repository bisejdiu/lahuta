#include <rdkit/GraphMol/RWMol.h>

#include "compute/context.hpp"
#include "compute/result.hpp"
#include "models/topology/data.hpp"
#include "models/topology/kernels.hpp"
#include "typing/types.hpp"

// clang-format off
namespace lahuta::models::topology {

template <typename DataT>
ComputationResult
ModelBuildKernel::execute(DataContext<DataT, Mut::ReadWrite> &context, const ModelBuildParams &params) {
  auto &data = context.data();

  try {
    canonicalize_bond_indices(data);

    data.mol = std::make_shared<RDKit::RWMol>(data.vertices, data.bonds, data.graph_type); // make mol

    if (!data.conf) throw std::runtime_error("No conformer available for building molecule.");
    if (data.mol->getNumAtoms() != data.conf->getNumAtoms() || !data.conf->is3D()) {
      throw std::runtime_error(data.mol->getNumAtoms() != data.conf->getNumAtoms()
          ? "Conformer atom count does not match molecule atom count"
          : "Conformer must be 3D before adding to molecule");
    }
    data.mol->addConformer(data.conf.release(), true);

    build_ringinfo(data);
    correct_disulfide_atom_properties(data);

    return ComputationResult(true);
  } catch (const std::exception &e) {
    return ComputationResult(ComputationError(std::string("Error in model atom typing: ") + e.what()));
  }
}

void ModelBuildKernel::canonicalize_bond_indices(ModelData &data) {
  for (size_t i = 0; i < data.bonds.size(); ++i) {
    if (data.bonds[i]) {
      data.bonds[i]->setIdx(static_cast<unsigned int>(i));
    }
  }
}

void ModelBuildKernel::build_ringinfo(ModelData &data) {
  if (!data.mol) return;

  data.aromatic_bond_indices.clear();
  data.aromatic_bond_indices.reserve(data.aromatic_atom_indices.size());

  for (const auto &ring_atoms : data.aromatic_atom_indices) {
    std::vector<int> ring_bonds;
    const size_t n = ring_atoms.size();
    if (n < 3) {
      data.aromatic_bond_indices.push_back(ring_bonds);
      continue;
    }

    for (size_t i = 0; i < n; ++i) {
      const int u = ring_atoms[i];
      const int v = ring_atoms[(i + 1) % n];
      const RDKit::Bond *b = data.mol->getBondBetweenAtoms(u, v);
      if (b) {
        ring_bonds.push_back(static_cast<int>(b->getIdx()));
      } else {
        ring_bonds.push_back(-1);
      }
    }
    data.aromatic_bond_indices.push_back(std::move(ring_bonds));
  }

  if (data.mol->getRingInfo()->isInitialized()) {
    data.mol->getRingInfo()->reset();
  }
  data.mol->getRingInfo()->initialize(RDKit::FIND_RING_TYPE_SYMM_SSSR);
  data.mol->getRingInfo()->addAllRings(data.aromatic_atom_indices, data.aromatic_bond_indices);
}

void ModelBuildKernel::correct_disulfide_atom_properties(ModelData &data) {
  if (!data.mol) return;

  for (const auto &pair : data.disulfide_pairs) {
    data.mol->getAtomWithIdx(pair.first )->setNumCompImplicitHs(0);
    data.mol->getAtomWithIdx(pair.second)->setNumCompImplicitHs(0);

    auto first_at  = static_cast<AtomType>(data.mol->getAtomWithIdx(pair.first )->getCompAtomType());
    auto second_at = static_cast<AtomType>(data.mol->getAtomWithIdx(pair.second)->getCompAtomType());
    data.mol->getAtomWithIdx(pair.first )->setCompAtomType(static_cast<int>(first_at  ^ AtomType::HbondDonor));
    data.mol->getAtomWithIdx(pair.second)->setCompAtomType(static_cast<int>(second_at ^ AtomType::HbondDonor));
  }
}

template ComputationResult ModelBuildKernel::execute<ModelData>(DataContext<ModelData, Mut::ReadWrite> &, const ModelBuildParams &);

} // namespace lahuta::models::topology
