#include "chemistry/group_typing.hpp"
#include "compute/context.hpp"
#include "compute/result.hpp"
#include "logging.hpp"
#include "topology/context.hpp"
#include "topology/kernels.hpp"

// clang-format off
namespace lahuta::topology {

template <typename DataT>
ComputationResult
SeedFromModelKernel::execute(DataContext<DataT, Mut::ReadWrite> &context, const SeedFromModelParams &params) {
  auto &data = context.data();

  try {
    data.residues->build();

    // CompAtomType is set by ModelTopology path
    data.atoms.clear();
    data.atoms.reserve(data.mol->getNumAtoms());
    for (const auto atom : data.mol->atoms()) {
      auto atom_type = static_cast<AtomType>(atom->getCompAtomType());
      data.atoms.push_back(AtomRec{atom_type, *atom});
    }

    // FIX: these are not fast. I'm working on alternatives.
    data.groups = GroupTypeAnalysis::analyze(*data.mol, *data.residues);
    data.rings  = AtomTypingKernel::populate_ring_entities(*data.mol);

    Logger::get_logger()->debug("seed_from_model: atoms={}, groups={}, rings={}", data.atoms.size(), data.groups.size(), data.rings.size());
    return ComputationResult(true);
  } catch (const std::exception &e) {
    return ComputationResult(ComputationError(std::string("Error seeding from model: ") + e.what()));
  }
}

template ComputationResult SeedFromModelKernel::execute<TopologyContext>(DataContext<TopologyContext, Mut::ReadWrite> &, const SeedFromModelParams &);

} // namespace lahuta::topology
