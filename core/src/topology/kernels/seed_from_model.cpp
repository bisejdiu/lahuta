#include "chemistry/group_typing.hpp"
#include "compute/context.hpp"
#include "compute/result.hpp"
#include "logging.hpp"
#include "topology/context.hpp"
#include "topology/kernels.hpp"
#include "typing/types.hpp"

// clang-format off
namespace lahuta::topology {

template <typename DataT>
ComputationResult
SeedFromModelKernel::execute(DataContext<DataT, Mut::ReadWrite> &context, const SeedFromModelParams &params) {
  auto &data = context.data();

  try {
    data.residues->build();

    data.atoms.clear();
    data.atoms.reserve(data.mol->getNumAtoms());
    if (params.mode == ContactComputerType::Molstar) {
      for (const auto atom : data.mol->atoms()) {
        auto atom_type = static_cast<AtomType>(atom->getCompAtomType());
        data.atoms.push_back(AtomRec{atom_type, *atom});
      }
    } else {
      for (const auto atom : data.mol->atoms()) {
        AtomType atom_type = get_atom_type(atom);
        data.atoms.push_back(AtomRec{atom_type, *atom});
      }
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
