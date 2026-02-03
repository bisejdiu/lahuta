/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto parts = std::make_tuple("besian", "sejdiu", "@gmail.com");
 *   return std::apply([](auto... p) { std::string s; (s.append(p), ...); return s; }, parts);
 * }();
 *
 */

#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>

#include "compute/context.hpp"
#include "compute/result.hpp"
#include "models/fast_lookup.hpp"
#include "models/tables.hpp"
#include "models/topology/data.hpp"
#include "models/topology/kernels.hpp"

// clang-format off
namespace lahuta::models::topology {
namespace C = lahuta::compute;

template <typename DataT>
C::ComputationResult
ModelAtomsKernel::execute(C::DataContext<DataT, C::Mut::ReadWrite> &context, const ModelAtomsParams &params) {
  auto &data = context.data();

  try {
    const std::string &seq = data.input_data->sequence;
    const size_t num_atoms = data.input_data->coords.size();
    const int num_residues = static_cast<int>(seq.size());

    data.graph_type = params.graph_type;

    data.vertices.clear();
    data.vertices.reserve(num_atoms);

    data.sulphur_atom_indices.clear();
    data.atom_idx = 0;

    for (int residue_idx = 0; residue_idx < num_residues; ++residue_idx) {
      const char aa_type = seq[residue_idx];
      create_atoms_for_residue(residue_idx, aa_type, data);
    }

    // terminal OXT atom
    if (!seq.empty()) {
      fix_termini_atoms(seq.back(), num_residues, data);
    }

    return C::ComputationResult(true);
  } catch (const std::exception &e) {
    return C::ComputationResult(C::ComputationError(std::string("Error creating model atoms: ") + e.what()));
  }
}

void ModelAtomsKernel::create_atoms_for_residue(int residue_idx, char aa_type, ModelData &data) {

  const auto &entry = StandardAminoAcidDataTable[aa_type];
  for (size_t local_atom_index = 0; local_atom_index < entry.size; ++local_atom_index) {
    const char *atom_name = entry.atoms[local_atom_index];
    const int   ih        = entry.ih[local_atom_index];
    const int   at        = entry.at[local_atom_index];
    const int   hyb       = entry.hyb[local_atom_index];
    const int   atom_num  = StandardAminoAcidAtomicNumbers[atom_name[0]];

    RDKit::Atom *atom = data.atom_pool->createAtom();
    atom->setIdx(static_cast<unsigned int>(data.atom_idx));
    atom->setAtomicNum(atom_num);
    double plddt_value = 0.0;
    const auto &plddts = data.input_data->plddt_per_residue;
    const auto res_idx = static_cast<std::size_t>(residue_idx);
    if (!plddts.empty() && res_idx < plddts.size()) {
      plddt_value = static_cast<double>(static_cast<std::uint8_t>(plddts[res_idx]));
    }
    atom->setMonomerInfo(data.info_pool->createAtomInfo(
        atom_name,
        data.atom_idx + 1,
        entry.name,
        residue_idx + 1,
        plddt_value));
    atom->setNumCompImplicitHs(ih);
    atom->setCompAtomType(at);
    atom->setHybridization(static_cast<RDKit::Atom::HybridizationType>(hyb));

    if (atom_num == 16 && aa_type == 'C') { // sulfur in CYS only
      data.sulphur_atom_indices.push_back(data.atom_idx);
    }

    data.vertices.push_back(atom);
    ++data.atom_idx;
  }
}

void ModelAtomsKernel::fix_termini_atoms(char aa_type, int num_residues, ModelData &data) {

  const auto entry = StandardAminoAcidDataTable[aa_type];

  // add OXT atom for the last residue
  RDKit::Atom *oxt_atom = data.atom_pool->createAtom();
  oxt_atom->setIdx(static_cast<unsigned int>(data.atom_idx));
  oxt_atom->setAtomicNum(8);
  double plddt_value = 0.0;
  const auto &plddts = data.input_data->plddt_per_residue;
  if (!plddts.empty() && num_residues > 0) {
    const auto idx = static_cast<std::size_t>(num_residues - 1);
    if (idx < plddts.size()) {
      plddt_value = static_cast<double>(static_cast<std::uint8_t>(plddts[idx]));
    }
  }
  oxt_atom->setMonomerInfo(data.info_pool->createAtomInfo(
      "OXT",
      data.atom_idx + 1,
      entry.name,
      num_residues,
      plddt_value));
  oxt_atom->setNumCompImplicitHs(0);
  oxt_atom->setCompAtomType(9217);
  oxt_atom->setHybridization(RDKit::Atom::SP2);

  data.vertices.push_back(oxt_atom);

  // first atom (N of first residue) gets SP3 and 3 implicit Hs
  if (!data.vertices.empty()) {
    data.vertices.front()->setHybridization(RDKit::Atom::SP3);
    data.vertices.front()->setNumCompImplicitHs(3);
  }

  ++data.atom_idx;
}

template C::ComputationResult ModelAtomsKernel::execute<ModelData>(C::DataContext<ModelData, C::Mut::ReadWrite> &, const ModelAtomsParams &);

} // namespace lahuta::models::topology
