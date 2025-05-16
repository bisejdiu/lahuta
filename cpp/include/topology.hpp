#ifndef LAHUTA_TOPOLOGY_HPP
#define LAHUTA_TOPOLOGY_HPP

#include "residues.hpp"
#include "topology/engine.hpp"
#include "topology_flags.hpp"
#include <memory>

namespace lahuta {

// FIX: using a "dynamic" cutoff might be better. For common atoms use a small cutoff. For other 
// atoms we'd use a larger cutoff but only around them.
constexpr static float BONDED_NEIGHBOR_SEARCH_CUTOFF = 4.5;
enum class ContactComputerType { None, Arpeggio, Molstar };

// Options for configuring topology parameters
struct TopologyBuildingOptions {
  ContactComputerType atom_typing_method = ContactComputerType::Molstar;
  double cutoff = BONDED_NEIGHBOR_SEARCH_CUTOFF;
  bool auto_heal = true; // auto-healing of dependencies
  bool compute_nonstandard_bonds = true; // whether to compute bonds for non-standard atoms
};

class Topology {
public:
  Topology() = default;
  Topology(std::shared_ptr<RDKit::RWMol> mol) 
    : mol_(mol), engine_(std::make_unique<topology::TopologyEngine>(mol)) {}

  const Residues &get_residues() const { return *engine_->get_data().residues; }
  const AtomEntityCollection  &get_atom_types() const { return engine_->get_data().atom_types; }
  const RingEntityCollection  &get_rings()      const { return engine_->get_data().rings; }
  const GroupEntityCollection &get_features()   const { return engine_->get_data().features; }

  std::vector<int> get_atom_ids() const { return get_residues().get_atom_ids(); }

  void build(TopologyBuildingOptions tops);

  void run_mask(TopologyComputation mask) const {
    for (auto bit : BASE_COMPUTATION_FLAGS)
      if (has_flag(mask, bit)) {
        engine_->get_engine()->run<void>(Topology::get_label(bit)); // auto-heal inside
      }
  }

  void assign_molstar_typing();
  void assign_arpeggio_atom_types();

  /// Enable/disable a specific computation
  void enable_computation(TopologyComputation comp, bool enabled);

  /// Enable only the specified computations (disabling all others)
  void enable_only(TopologyComputation comps);

  /// Check if a specific computation is enabled
  bool is_computation_enabled(TopologyComputation comp) const;

  /// Execute a specific computation (with dependencies)
  bool execute_computation(TopologyComputation comp);

  /// Set the neighbor search cutoff
  void set_cutoff(double cutoff);

  /// Set the atom typing method
  void set_atom_typing_method(ContactComputerType method);

  /// Set whether to compute non-standard bonds
  void set_compute_nonstandard_bonds(bool compute);

  /// Get the engine
  topology::TopologyEngine* get_engine() { return engine_.get(); }

  /// approximate total memory usage
  size_t total_size() const;

private:
  // Get compute::ComputationLabel from TopologyComputation
  static const topology::ComputationLabel& get_label(TopologyComputation comp);

private:
  std::shared_ptr<RDKit::RWMol> mol_;
  std::unique_ptr<topology::TopologyEngine> engine_;
};

} // namespace lahuta

#endif // LAHUTA_TOPOLOGY_HPP
