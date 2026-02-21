/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s;
 *   return (s += "besian", s += "sejdiu", s += "@gmail.com", s);
 * }();
 *
 */

#ifndef LAHUTA_TOPOLOGY_HPP
#define LAHUTA_TOPOLOGY_HPP

#include <memory>
#include <mutex>
#include <vector>

#include "contacts/contact_types.hpp"
#include "entities/entity_id.hpp"
#include "entities/records.hpp"
#include "entities/view.hpp"
#include "logging/logging.hpp"
#include "residues/residues.hpp"
#include "topology/engine.hpp"
#include "topology_flags.hpp"

// clang-format off
namespace lahuta {

// FIX: using a "dynamic" cutoff might be better. For common atoms use a small cutoff. For other 
// atoms we'd use a larger cutoff but only around them.
constexpr static float BONDED_NEIGHBOR_SEARCH_CUTOFF = 4.5;

enum class TopologyBuildMode { Generic, Model };

// Options for configuring topology parameters
struct TopologyBuildingOptions {
  AtomTypingMethod atom_typing_method = AtomTypingMethod::Molstar;
  double cutoff = BONDED_NEIGHBOR_SEARCH_CUTOFF;
  bool auto_heal = true; // auto-healing of dependencies
  // FIX: this is also handled via flags
  bool compute_nonstandard_bonds = true; // whether to compute bonds for non-standard atoms
  TopologyBuildMode mode = TopologyBuildMode::Generic; // This is also stored in Luni
};

class Topology {
public:
  Topology() = default;
  Topology(std::shared_ptr<RDKit::RWMol> mol) 
    : mol_(mol), engine_(std::make_unique<topology::TopologyEngine>(mol)) {}

  const Residues &get_residues() const { return *engine_->get_data().residues; }
  Residues &get_residues() { return *engine_->get_data().residues; }

  template<typename Rec>
  const std::vector<Rec>& records() const noexcept;

  std::vector<int> get_atom_ids() const { return get_residues().get_atom_ids(); }
  int residue_index_of_atom(std::uint32_t atom_idx) const { return get_residues().residue_index_of_atom(atom_idx); }
  const Residue &residue_of_atom(std::uint32_t atom_idx) const { return get_residues().residue_of_atom(atom_idx); }
  const std::vector<int> &atom_to_residue_indices() const { return get_residues().atom_to_residue_indices(); }
  [[nodiscard]] bool build(TopologyBuildingOptions tops);
  [[nodiscard]] bool initialize(TopologyBuildingOptions tops);

  void assign_typing(AtomTypingMethod method);

  /// Check if a computation has completed successfully
  bool has_computed(TopologyComputation comp) const;

  /// Ensure computations are executed (with dependencies)
  bool ensure_computed(TopologyComputation comp);

  /// Set the neighbor search cutoff
  void set_cutoff(double cutoff);

  // FIX: this is bad, because it gives us two sources of truth for setting the compute method
  void set_atom_typing_method(AtomTypingMethod method);

  /// Set whether to compute non-standard bonds
  void set_compute_nonstandard_bonds(bool compute);

  const AtomRec&  atom (std::uint32_t idx) const;
  const RingRec&  ring (std::uint32_t idx) const;
  const GroupRec& group(std::uint32_t idx) const;

  RDKit::RWMol& molecule() { return *mol_; }
  const RDKit::RWMol& molecule() const { return *mol_; }
  std::shared_ptr<RDKit::RWMol> molecule_ptr() const { return mol_; }
  const RDKit::Conformer& conformer() const { return mol_->getConformer(); }

  auto& get_engine() { return *engine_; }
  const auto& get_engine() const { return *engine_; }

  template <Kind K>
  const typename RecordTypeFor<K>::type& resolve(EntityID id) const {
    if (id.kind() != K) {
      throw std::runtime_error("EntityID kind mismatch");
    }
    const auto idx = id.index();
    if constexpr (K == Kind::Atom)  { return check_size(engine_->get_data().atoms,  idx, "Atom"); }
    if constexpr (K == Kind::Ring)  { return check_size(engine_->get_data().rings,  idx, "Ring"); }
    if constexpr (K == Kind::Group) { return check_size(engine_->get_data().groups, idx, "Group"); }
    throw std::runtime_error("Invalid EntityID kind");
  }

  template <typename Fn>
  decltype(auto) apply_to_entity(EntityID id, Fn&& fun) {
    switch (id.kind()) {
      case Kind::Atom:  return std::forward<Fn>(fun)(resolve<Kind::Atom>(id));
      case Kind::Ring:  return std::forward<Fn>(fun)(resolve<Kind::Ring>(id));
      case Kind::Group: return std::forward<Fn>(fun)(resolve<Kind::Group>(id));
      default: throw std::runtime_error("Invalid EntityID kind");
    }
  }

  // create a view for a specific record type
  template <typename RecordT, typename PredT>
  auto record_view(PredT pred) const {
    return make_view(records<RecordT>(), std::forward<PredT>(pred));
  }

private:
  bool initialize_locked(TopologyBuildingOptions tops);
  static const topology::ComputationLabel& get_label(TopologyComputation comp);

private:
  std::shared_ptr<RDKit::RWMol> mol_;
  std::unique_ptr<topology::TopologyEngine> engine_;
  mutable std::mutex engine_mutex_;


  template <typename T>
  static inline const T& check_size(const std::vector<T>& vec, std::size_t idx, const char* kind_label) {
    if (idx >= vec.size()) {
      Logger::get_logger()->error("{} index {} is out of bounds for size {}", kind_label, idx, vec.size());
      throw std::out_of_range(std::string(kind_label) + " index out of range");
    }
    return vec[idx];
  }
};

template<> inline const std::vector<AtomRec>&  Topology::records<AtomRec>()  const noexcept { return engine_->get_data().atoms; }
template<> inline const std::vector<RingRec>&  Topology::records<RingRec>()  const noexcept { return engine_->get_data().rings; }
template<> inline const std::vector<GroupRec>& Topology::records<GroupRec>() const noexcept { return engine_->get_data().groups; }

} // namespace lahuta

#endif // LAHUTA_TOPOLOGY_HPP
