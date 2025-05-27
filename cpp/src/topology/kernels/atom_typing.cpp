#include "compute/context.hpp"
#include "compute/result.hpp"
#include "contacts/atoms.hpp"
#include "contacts/groups.hpp"
#include "logging.hpp"
#include "topology/data.hpp"
#include "topology/kernels.hpp"
#include <rdkit/GraphMol/RWMol.h>
#include <valence_model.hpp>

// clang-format off
namespace lahuta::topology {

template <typename DataT>
ComputationResult
AtomTypingKernel::execute(DataContext<DataT, Mut::ReadWrite> &context, const AtomTypingParams &params) {
  auto &data = context.data();

  try {
    ValenceModel valence_model;
    valence_model.apply(*data.mol);

    data.atoms  = AtomTypeAnalysis()(*data.mol);
    data.groups = GroupTypeAnalysis::analyze(*data.mol, *data.residues);
    data.rings  = populate_ring_entities(*data.mol);

    return ComputationResult(true);
  } catch (const std::exception &e) {
    Logger::get_logger()->error("Exception in atom typing: {}", e.what());
    return ComputationResult(ComputationError(std::string("Error computing atom types: ") + e.what()));
  }
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
    std::vector<std::uint32_t> atom_indices;
    atom_indices.reserve(ring.size());

    for (int atom_idx : ring) {
      atom_indices.push_back(static_cast<std::uint32_t>(atom_idx));
    }

    // center of the ring
    RDGeom::Point3D center{0.0, 0.0, 0.0};
    for (int atom_idx : ring) {
      center += conf.getAtomPos(atom_idx);
    }
    center /= ring.size();

    // FIX: bad impl normal vector
    RDGeom::Point3D normal{0.0, 0.0, 1.0};
    if (ring.size() >= 3) {
      auto p1 = conf.getAtomPos(ring[0]);
      auto p2 = conf.getAtomPos(ring[1]);
      auto p3 = conf.getAtomPos(ring[2]);

      auto v1 = p2 - p1;
      auto v2 = p3 - p1;
      normal = v1.crossProduct(v2);
      normal.normalize();
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
      /*.atoms    =*/ std::move(atom_indices),
      /*.center   =*/ center,
      /*.normal   =*/ normal,
      /*.aromatic =*/ is_aromatic
    });
  }

  return ring_recs;
}

bool AtomTypingKernel::should_initialize_ringinfo(int mol_size) {
  constexpr int small_threshold  = 20'000;
  constexpr int medium_threshold = 50'000;
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

template ComputationResult AtomTypingKernel::execute<TopologyData>(DataContext<TopologyData, Mut::ReadWrite> &, const AtomTypingParams &);

} // namespace lahuta::topology
