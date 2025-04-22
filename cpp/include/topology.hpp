#ifndef LAHUTA_TOPOLOGY_HPP
#define LAHUTA_TOPOLOGY_HPP

#include "contacts/atoms.hpp"
#include "contacts/groups.hpp"
#include "residues.hpp"

// clang-format off
namespace lahuta {

// FIX: using a "dynamic" cutoff might be better. For common atoms use a small cutoff. For other 
// atoms we'd use a larger cutoff but only around them.
constexpr static float BONDED_NEIGHBOR_SEARCH_CUTOFF = 4.5;
enum class ContactComputerType { None, Arpeggio, Molstar };

struct TopologyBuildingOptions {
  const bool identify_ring_atoms = true; // if we decide to also check for ring atoms using atom names:e.g. only for protein-only systems as an optimization technique
  ContactComputerType atom_typing_method = ContactComputerType::Molstar;
  double cutoff = BONDED_NEIGHBOR_SEARCH_CUTOFF;
  bool compute_bonds = true;
};


class Topology {
public:
  Topology() = default;
  Topology(std::shared_ptr<RDKit::RWMol> mol) : mol_(mol), residues(std::make_unique<Residues>(*mol)) {}

  const Residues &get_residues() const { return *residues; }
  const AtomEntityCollection  &get_atom_types() const { return atom_types; }
  const RingEntityCollection  &get_rings()      const { return rings_vec; }
  const GroupEntityCollection &get_features()   const { return features; }

  std::vector<int> get_atom_ids() const { return residues->get_atom_ids(); }

  void build(TopologyBuildingOptions tops);

  void assign_molstar_typing() {

    ValenceModel valence_model;
    valence_model.apply(*mol_);

    // FIX: to be replaced by the entitytype manager
    atom_types = AtomTypeAnalysis ::analyze(*mol_);
    features   = GroupTypeAnalysis::analyze(*mol_, *residues);

    // FIX: We have a conceptual issue with Aromaticity:
    //  - `features` includes AromaticRingGroups.
    //  - `rings_vec` includes ***only*** aromatic rings
    //  - `atom_types`, which stores atom-level identifiers, does not store their aromaticity identifier

    rings_vec = populate_ring_entities();
  }

  void assign_arpeggio_atom_types();

  /// approximate total memory usage
  size_t total_size() const;

private:
  RingEntityCollection populate_ring_entities();
  void compute_bonds(const NSResults &neighbors);
  static void cleanup_predef(RDKit::RWMol &mol);
  static void cleanup(RDKit::RWMol &mol);
  static void merge_bonds(RDKit::RWMol &target, RDKit::RWMol &source, const std::vector<int> &index_map);
  static bool should_initialize_ringinfo(int mol_size);

private:
  AtomEntityCollection  atom_types;
  RingEntityCollection  rings_vec;
  GroupEntityCollection features;

  std::shared_ptr<RDKit::RWMol> mol_;
  std::unique_ptr<Residues> residues;
};

} // namespace lahuta

#endif // LAHUTA_TOPOLOGY_HPP
