#include <models/tables.hpp>
#include <rdkit/GraphMol/RWMol.h>

#include "compute/context.hpp"
#include "compute/result.hpp"
#include "models/constants.hpp"
#include "models/topology/data.hpp"
#include "models/topology/kernels.hpp"

// clang-format off
namespace lahuta::models::topology {
namespace C = lahuta::compute;

template <typename DataT>
C::ComputationResult
ModelAromaticsKernel::execute(C::DataContext<DataT, C::Mut::ReadWrite> &context, const ModelAromaticsParams &params) {
  auto &data = context.data();

  try {
    const std::string &seq = data.input_data->sequence;
    const int num_residues = static_cast<int>(seq.size());

    data.aromatic_atom_indices.clear();
    data.aromatic_bond_indices.clear();

    int atom_idx = 0;
    for (int residue_idx = 0; residue_idx < num_residues; ++residue_idx) {
      const char aa_type = seq[residue_idx];
      const auto &entry  = StandardAminoAcidDataTable[aa_type];
      const int residue_start_idx = atom_idx;

      switch (aa_type) {
        case 'F': process_aromatic_residue(residue_start_idx, phe_arom_indices,  data); break;
        case 'Y': process_aromatic_residue(residue_start_idx, tyr_arom_indices,  data); break;
        case 'H': process_aromatic_residue(residue_start_idx, his_arom_indices,  data); break;
        case 'W':
          process_aromatic_residue(residue_start_idx, trp_arom_indices5, data);
          process_aromatic_residue(residue_start_idx, trp_arom_indices6, data);
          break;
        default:
          break;
      }

      atom_idx += static_cast<int>(entry.size);
    }

    return C::ComputationResult(true);
  } catch (const std::exception &e) {
    return C::ComputationResult(C::ComputationError(std::string("Error processing model aromatics: ") + e.what()));
  }
}

template <typename ArrayN>
void ModelAromaticsKernel::process_aromatic_residue(int residue_start_idx, const ArrayN &atom_indices, ModelData &data) {

  std::vector<int> ring_atom_indices;
  ring_atom_indices.reserve(atom_indices.size());

  for (int idx : atom_indices) {
    const int global_idx = residue_start_idx + idx;
    ring_atom_indices.push_back(global_idx);
    if (global_idx < static_cast<int>(data.vertices.size())) {
      data.vertices[static_cast<size_t>(global_idx)]->setIsAromatic(true);
    }
  }

  data.aromatic_atom_indices.push_back(std::move(ring_atom_indices));
}

template C::ComputationResult ModelAromaticsKernel::execute<ModelData>(C::DataContext<ModelData, C::Mut::ReadWrite> &, const ModelAromaticsParams &);

} // namespace lahuta::models::topology
