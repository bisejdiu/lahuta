#ifndef LAHUTA_TOPOLOGY_HPP
#define LAHUTA_TOPOLOGY_HPP

#include "residues.hpp"
#include "topology/engine.hpp"
#include "topology_flags.hpp"
#include "entities/entity_id.hpp"
#include "entities/records.hpp"
#include "entities/view.hpp"
#include <memory>
#include <vector>

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
  Residues &get_residues() { return *engine_->get_data().residues; }

  template<typename Rec>
  constexpr const std::vector<Rec>& records() const noexcept;

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

  // FIX: this is bad, because it gives us two sources of truth for setting the compute method
  void set_atom_typing_method(ContactComputerType method);

  /// Set whether to compute non-standard bonds
  void set_compute_nonstandard_bonds(bool compute);

  /// approximate total memory usage
  size_t total_size() const;

  const AtomRec&  atom (std::uint32_t idx) const;
  const RingRec&  ring (std::uint32_t idx) const;
  const GroupRec& group(std::uint32_t idx) const;

  RDKit::RWMol& molecule() { return *mol_; }
  const RDKit::RWMol& molecule() const { return *mol_; }
  const RDKit::Conformer& conformer() const { return mol_->getConformer(); }

  auto& get_engine() { return *engine_; }
  const auto& get_engine() const { return *engine_; }

  template <Kind K>
  const typename RecordTypeFor<K>::type& resolve(EntityID id) const {
    if (id.kind() != K) {
      throw std::runtime_error("EntityID kind mismatch");
    }
    switch (K) {
      case Kind::Atom:  return atom(id.index());
      case Kind::Ring:  return ring(id.index());
      case Kind::Group: return group(id.index());
      default: throw std::runtime_error("Invalid EntityID kind");
    }
  }

  template <typename Fn>
  decltype(auto) apply_to_entity(EntityID id, Fn&& fun) {
    switch (id.kind()) {
      case Kind::Atom:  return fun(atom (id.index()));
      case Kind::Ring:  return fun(ring (id.index()));
      case Kind::Group: return fun(group(id.index()));
      default: throw std::runtime_error("Invalid EntityID kind");
    }
  }

  // create a view for a specific record type
  template <typename RecordT, typename PredT>
  auto record_view(PredT pred) const {
    return make_view(records<RecordT>(), std::forward<PredT>(pred));
  }

private:
  static const topology::ComputationLabel& get_label(TopologyComputation comp);

private:
  std::shared_ptr<RDKit::RWMol> mol_;
  std::unique_ptr<topology::TopologyEngine> engine_;
};

template<> inline constexpr const std::vector<AtomRec>&  Topology::records<AtomRec>()  const noexcept { return engine_->get_data().atoms; }
template<> inline constexpr const std::vector<RingRec>&  Topology::records<RingRec>()  const noexcept { return engine_->get_data().rings; }
template<> inline constexpr const std::vector<GroupRec>& Topology::records<GroupRec>() const noexcept { return engine_->get_data().groups; }

} // namespace lahuta

#endif // LAHUTA_TOPOLOGY_HPP
