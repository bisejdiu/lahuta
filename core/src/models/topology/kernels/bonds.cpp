/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, '\0'); std::size_t pos = 0;
 *   auto add = [&](std::string_view p) { std::copy(p.begin(), p.end(), &s[pos]); pos += p.size(); };
 *   add("besian"); add("sejdiu"); add("@"); add("gmail.com");
 *   return s;
 * }();
 *
 */

#include <rdkit/GraphMol/RWMol.h>

#include "compute/context.hpp"
#include "compute/result.hpp"
#include "models/tables.hpp"
#include "models/topology/data.hpp"
#include "models/topology/kernels.hpp"

namespace lahuta::models::topology {
namespace C = lahuta::compute;

template <typename DataT>
C::ComputationResult ModelBondsKernel::execute(C::DataContext<DataT, C::Mut::ReadWrite> &context,
                                               const ModelBondsParams &params) {
  auto &data = context.data();

  try {
    const std::string &seq = data.input_data->sequence;
    const int num_residues = static_cast<int>(seq.size());

    // estimate bond count
    size_t expected_bond_count = 0;
    for (int residue_idx = 0; residue_idx < num_residues; ++residue_idx) {
      const char aa_type           = seq[residue_idx];
      const AminoAcidEdges &edges  = StandardAminoAcidBondTable[aa_type];
      expected_bond_count         += edges.size;
    }

    data.bonds.clear();
    data.bonds.reserve(expected_bond_count);

    int atom_idx = 0;
    for (int residue_idx = 0; residue_idx < num_residues; ++residue_idx) {
      const char aa_type          = seq[residue_idx];
      const auto &entry           = StandardAminoAcidDataTable[aa_type];
      const int residue_start_idx = atom_idx;

      // intra-residue bonds
      const AminoAcidEdges &edges = StandardAminoAcidBondTable[aa_type];
      add_intra_residue_bonds(edges, residue_start_idx, data);

      // peptide bond to next residue
      atom_idx += static_cast<int>(entry.size);
      add_peptide_bond(residue_idx, num_residues, residue_start_idx, atom_idx, data);
    }

    // C-OXT bond
    if (!seq.empty()) {
      add_terminal_oxt_bond(seq.back(), data);
    }

    return C::ComputationResult(true);
  } catch (const std::exception &e) {
    return C::ComputationResult(C::ComputationError(std::string("Error creating model bonds: ") + e.what()));
  }
}

void ModelBondsKernel::add_intra_residue_bonds(const AminoAcidEdges &edges, int residue_start_idx,
                                               ModelData &data) {

  for (size_t j = 0; j < edges.size; ++j) {
    const auto &edge    = edges.edges[j];
    const int atom1_idx = residue_start_idx + edge.i;
    const int atom2_idx = residue_start_idx + edge.j;
    RDKit::Bond *bond   = data.bond_pool->createBond(atom1_idx, atom2_idx, edge.order);
    data.bonds.push_back(bond);
  }
}

void ModelBondsKernel::add_peptide_bond(int residue_idx, int num_residues, int residue_start_idx,
                                        int next_residue_start_idx, ModelData &data) {

  if (residue_idx < num_residues - 1) {
    const int current_c_idx = residue_start_idx + 2;  // C is always at position 2
    const int next_n_idx    = next_residue_start_idx; // start index of next residue (N at position 0)
    RDKit::Bond *bond       = data.bond_pool->createBond(current_c_idx, next_n_idx, RDKit::Bond::SINGLE);
    data.bonds.push_back(bond);
  }
}

void ModelBondsKernel::add_terminal_oxt_bond(char aa_type, ModelData &data) {

  const auto entry = StandardAminoAcidDataTable[aa_type];
  // C is at position 2, -1 because atom_idx was incremented
  const int last_residue_c_idx = data.atom_idx - static_cast<int>(entry.size) + 2 - 1;
  const int oxt_idx            = data.atom_idx - 1; // OXT is the last atom added
  RDKit::Bond *c_oxt_bond      = data.bond_pool->createBond(last_residue_c_idx, oxt_idx, RDKit::Bond::SINGLE);
  data.bonds.push_back(c_oxt_bond);
}

template C::ComputationResult
ModelBondsKernel::execute<ModelData>(C::DataContext<ModelData, C::Mut::ReadWrite> &,
                                     const ModelBondsParams &);

} // namespace lahuta::models::topology
