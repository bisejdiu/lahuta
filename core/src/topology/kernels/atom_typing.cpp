#include <algorithm>
#include <functional>
#include <optional>

#include <rdkit/GraphMol/RWMol.h>

#include "chemistry/atom_typing.hpp"
#include "chemistry/group_typing.hpp"
#include "chemistry/valence_model.hpp"
#include "compute/context.hpp"
#include "compute/result.hpp"
#include "logging/logging.hpp"
#include "residues/residues.hpp"
#include "selections/mol_filters.hpp"
#include "topology/context.hpp"
#include "topology/kernels.hpp"
#include "typing/flags.hpp"
#include "typing/getcontacts/atom_typing.hpp"
#include "typing/smarts_matching.hpp"
#include "typing/types.hpp"

namespace {

class ResidueChunkIterator {
public:
  ResidueChunkIterator(const lahuta::Residues &residues, int chunk_threshold)
      : residues_(residues.get_residues()), threshold_(std::max(chunk_threshold, 1)) {}

  std::optional<std::vector<int>> next() {
    while (next_residue_index_ < residues_.size()) {
      const auto &residue = residues_[next_residue_index_++];
      if (residue.atoms.empty()) continue;

      buffer_.reserve(buffer_.size() + residue.atoms.size());
      for (const auto *atom : residue.atoms) {
        buffer_.push_back(atom->getIdx());
      }

      if (static_cast<int>(buffer_.size()) >= threshold_) {
        return flush_buffer();
      }
    }

    if (buffer_.empty()) return std::nullopt;
    return flush_buffer();
  }

private:
  std::optional<std::vector<int>> flush_buffer() {
    std::vector<int> chunk;
    chunk.swap(buffer_);
    return chunk;
  }

  const std::vector<lahuta::Residue> &residues_;
  size_t next_residue_index_ = 0;
  int threshold_;
  std::vector<int> buffer_;
};

} // namespace

// clang-format off
namespace lahuta::topology {

template <typename DataT>
ComputationResult
AtomTypingKernel::execute(DataContext<DataT, Mut::ReadWrite> &context, const AtomTypingParams &params) {
  auto &data = context.data();

  try {
    switch (params.mode) {
      case AtomTypingMethod::Molstar: {
        ValenceModel valence_model;
        valence_model.apply(*data.mol);

        data.atoms  = AtomTypeAnalysis()(*data.mol);
        data.groups = GroupTypeAnalysis::analyze(*data.mol, *data.residues);

        Logger::get_logger()->debug("atom_typing: molstar, atoms={}, groups={}, rings={}", data.atoms.size(), data.groups.size(), data.rings.size());
        return ComputationResult(true);
      }
      case AtomTypingMethod::Arpeggio: {
        data.atoms.clear(); data.groups.clear();
        data.atoms.reserve(data.mol->getNumAtoms());
        for (auto atom : data.mol->atoms()) {
          AtomType atom_type = get_atom_type(atom);
          data.atoms.push_back(AtomRec{
            /*.type =*/  atom_type,
            /*,.atom =*/ *atom
          });
        }

        constexpr int residue_chunk_soft_threshold = 5'000;
        auto unknown_residues = data.residues->filter(std::not_fn(definitions::is_protein_extended));
        ResidueChunkIterator chunker(unknown_residues, residue_chunk_soft_threshold);

        while (auto chunk = chunker.next()) {
          if (chunk->empty()) continue;

          auto chunk_mol = filter_with_bonds(*data.mol, *chunk);
          if (!should_initialize_ringinfo(chunk_mol.getNumAtoms())) continue;

          chunk_mol.getRingInfo()->initialize(RDKit::FIND_RING_TYPE_SYMM_SSSR);
          RDKit::MolOps::findSSSR(chunk_mol);

          auto chunk_types = match_atom_types(chunk_mol);
          for (size_t i = 0; i < chunk->size(); ++i) {
            auto &rec = data.atoms[(*chunk)[i]];
            rec.type |= chunk_types[i];
          }
        }

        // Set aromatic flags on atoms based on rings computed by RingComputation dependency
        for (const auto &ring : data.rings) {
          if (ring.aromatic) {
            for (const auto &atom_ref : ring.atoms) {
              auto &atom = atom_ref.get();
              auto &rec = data.atoms[atom.getIdx()];
              rec.type |= AtomType::Aromatic;
            }
          }
        }
        Logger::get_logger()->debug("atom_typing: arpeggio, atoms={}, groups={}, rings={}", data.atoms.size(), data.groups.size(), data.rings.size());
        return ComputationResult(true);
      }
      case AtomTypingMethod::GetContacts: {
        data.atoms.clear(); data.groups.clear();
        data.atoms.reserve(data.mol->getNumAtoms());
        for (auto atom : data.mol->atoms()) {
          AtomType atom_type = typing::getcontacts::classify_atom(*data.mol, *atom);
          data.atoms.push_back(AtomRec{
            /*.type =*/  atom_type,
            /*,.atom =*/ *atom
          });
        }

        // Set aromatic flags on atoms based on rings computed by RingComputation dependency
        for (const auto &ring : data.rings) {
          if (ring.aromatic) {
            for (const auto &atom_ref : ring.atoms) {
              auto &atom = atom_ref.get();
              data.atoms[atom.getIdx()].type |= AtomType::Aromatic;
            }
          }
        }

        Logger::get_logger()->debug("atom_typing: getcontacts, atoms={}, groups={}, rings={}", data.atoms.size(), data.groups.size(), data.rings.size());
        return ComputationResult(true);
      }
    }
  } catch (const std::exception &e) {
    Logger::get_logger()->error("Exception in atom typing: {}", e.what());
    return ComputationResult(ComputationError(std::string("Error computing atom types: ") + e.what()));
  }
  return ComputationResult(ComputationError("Unknown atom typing method"));
}

//
// Ring perception and RingRec creation are unnaturally separated. Our choice to use RDKit's RingInfo,
// which was done for efficient querying of ring information, forces us to treat rings differently. A further
// complication is that our GroupTypeAnalysis also stores rings!  - besian, May, 22, 2025
//
std::vector<RingRec> AtomTypingKernel::populate_ring_entities(RDKit::RWMol &mol) {
  std::vector<RingRec> ring_recs;

  const auto& conf  = mol.getConformer();
  const auto& rings = mol.getRingInfo()->atomRings();
  ring_recs.reserve(rings.size());

  for (const std::vector<int> &ring : rings) {
    std::vector<std::reference_wrapper<const RDKit::Atom>> atoms;
    std::vector<std::uint32_t> atom_indices;
    atom_indices.reserve(ring.size());
    atoms.reserve(ring.size());

    for (int atom_idx : ring) {
      atom_indices.push_back(static_cast<std::uint32_t>(atom_idx));
      atoms.push_back(std::cref(*mol.getAtomWithIdx(atom_idx)));
    }

    // aromaticity check
    bool is_aromatic = false;
    for (int idx : ring) {
      if (mol.getAtomWithIdx(idx)->getIsAromatic()) {
        is_aromatic = true;
        break;
      }
    }

    ring_recs.push_back(RingRec{
      /*.atoms    =*/ std::move(atoms),
      /*.aromatic =*/ is_aromatic
    });
  }

  return ring_recs;
}

bool AtomTypingKernel::should_initialize_ringinfo(int mol_size) {
  constexpr int small_threshold  = 30'000;
  constexpr int medium_threshold = 65'000;
  constexpr int large_threshold  = 100'000;

  if (mol_size < small_threshold) return true;

  if (mol_size < medium_threshold) {
    Logger::get_logger()->warn("Filtered molecule size ({}) is large. Performance may be affected.", mol_size);
    return true;
  } else if (mol_size < large_threshold) {
    Logger::get_logger()->warn("Filtered molecule size ({}) is very large. Performance may be severely affected.", mol_size);
    return true;
  }

  Logger::get_logger()->error("Filtered molecule size ({}) is too large. Ring perception will be skipped!", mol_size);
  return false;
}

template ComputationResult AtomTypingKernel::execute<TopologyContext>(DataContext<TopologyContext, Mut::ReadWrite> &, const AtomTypingParams &);

} // namespace lahuta::topology
