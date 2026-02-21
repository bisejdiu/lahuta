/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto forward_concat = [](auto&& a, auto&& b, auto&& c) {
 *     return std::string(std::forward<decltype(a)>(a)) +
 *            std::forward<decltype(b)>(b) +
 *            std::forward<decltype(c)>(c);
 *   };
 *   return forward_concat("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_MODELS_TOPOLOGY_KERNELS_HPP
#define LAHUTA_MODELS_TOPOLOGY_KERNELS_HPP

#include <rdkit/GraphMol/RWMol.h>

#include "compute/context.hpp"
#include "compute/result.hpp"
#include "models/topology/data.hpp"
#include "models/topology/parameters.hpp"

namespace lahuta {
struct AminoAcidEdges;
} // namespace lahuta

// clang-format off
namespace lahuta::models::topology {
namespace C = lahuta::compute;

struct ModelAtomsKernel {
  template <typename DataT>
  static C::ComputationResult
  execute(C::DataContext<DataT, C::Mut::ReadWrite> &context, const ModelAtomsParams &params);

private:
  static void create_atoms_for_residue(int residue_idx, char aa_type, ModelData &data);
  static void fix_termini_atoms(char aa_type, int num_residues, ModelData &data);
};

struct ModelBondsKernel {
  template <typename DataT>
  static C::ComputationResult
  execute(C::DataContext<DataT, C::Mut::ReadWrite> &context, const ModelBondsParams &params);

private:
  static void add_intra_residue_bonds(const class AminoAcidEdges &edges, int residue_start_idx, ModelData &data);
  static void add_peptide_bond(int residue_idx, int num_residues, int residue_start_idx, int next_residue_start_idx, ModelData &data);
  static void add_terminal_oxt_bond(char aa_type, ModelData &data);
};

struct ModelPositionsKernel {
  template <typename DataT>
  static C::ComputationResult
  execute(C::DataContext<DataT, C::Mut::ReadWrite> &context, const ModelPositionsParams &params);
};

struct ModelAromaticsKernel {
  template <typename DataT>
  static C::ComputationResult
  execute(C::DataContext<DataT, C::Mut::ReadWrite> &context, const ModelAromaticsParams &params);

private:
  template <typename ArrayN>
  static void process_aromatic_residue(int residue_start_idx, const ArrayN &atom_indices, ModelData &data);
};

struct ModelDisulfidesKernel {
  template <typename DataT>
  static C::ComputationResult
  execute(C::DataContext<DataT, C::Mut::ReadWrite> &context, const ModelDisulfidesParams &params);
};

struct ModelBuildKernel {
  template <typename DataT>
  static C::ComputationResult
  execute(C::DataContext<DataT, C::Mut::ReadWrite> &context, const ModelBuildParams &params);

private:
  static void canonicalize_bond_indices(ModelData &data);
  static void build_ringinfo(ModelData &data);
  static void correct_disulfide_atom_properties(ModelData &data);
};

} // namespace lahuta::models::topology

#endif // LAHUTA_MODELS_TOPOLOGY_KERNELS_HPP
