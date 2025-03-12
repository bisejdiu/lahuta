#ifndef LAHUTA_TOPOLOGY_HPP
#define LAHUTA_TOPOLOGY_HPP

#include "aromatics.hpp"
#include "atom_types.hpp"
#include "contacts/atoms.hpp"
#include "contacts/groups.hpp"
#include "convert.hpp"
#include "definitions.hpp"
#include "residues.hpp"
#include "spdlog/spdlog.h"
#include <rdkit/GraphMol/BondIterators.h>

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


  void build(TopologyBuildingOptions tops) {

    if (!mol_) {
      spdlog::critical("Cannot build topology without a molecule.");
      throw std::runtime_error("Make sure to provide a molecule before building the topology.");
    }

    try {
      if (tops.compute_bonds) {
      auto grid = FastNS(mol_->getConformer().getPositions());
      auto ok = grid.build(tops.cutoff);
      if (!ok) {
          throw std::runtime_error("Failed to build the grid for neighbor search.");
      }
      auto neighbors = std::make_shared<NSResults>(grid.self_search());

      // FIX: neighbor computation can technically be the responsibility of compute_bonds, but
      // that moves them perhaps too much down the stack, and makes control of the process
      // more difficult (e.g., if we need to use a memory pool or arena allocator)
      this->compute_bonds(*neighbors);
      }

      // build residue information
      residues->build();

      // populate ring information to RDKit Mol
      initialize_and_populate_ringinfo(*mol_, *residues);

      switch (tops.atom_typing_method) {
        case ContactComputerType::Molstar:
          this->assign_molstar_typing();
          break;
        case ContactComputerType::Arpeggio:
          this->assign_arpeggio_atom_types();
          break;
        default:
          break;
      }

    } catch (const std::runtime_error &e) {
      spdlog::critical("Error creating topology! Exception caught: {}. Will not terminate, "
                       "but no topology-based features will be available.", e.what());
    }
  }

  void assign_molstar_typing() {
    std::cout << "Assigning MolStar atom types" << std::endl;

    // FIX: to be replaced by the entitytype manager
    atom_types = AtomTypeAnalysis ::analyze(*mol_);
    features   = GroupTypeAnalysis::analyze(*mol_, *residues);

    // FIX: We have a conceptual issue with Aromaticity:
    //  - `features` includes AromaticRingGroups.
    //  - `rings_vec` includes ***only*** aromatic rings
    //  - `atom_types`, which stores atom-level identifiers, does not store their aromaticity identifier

    rings_vec = populate_ring_entities();
  }

  void assign_arpeggio_atom_types() {

    std::cout << "Assigning Arpeggio atom types" << std::endl;

    atom_types.reserve(mol_->getNumAtoms());
    for (auto atom : mol_->atoms()) {
      AtomType atom_type = get_atom_type(atom);
      atom_types.add_data(*mol_, atom, atom_type);
    }

    auto unk_indices = residues->filter(std::not_fn(definitions::is_protein_extended)).get_atom_ids();
    if (!unk_indices.empty()) {
      std::sort(unk_indices.begin(), unk_indices.end());
      auto new_mol = filter_with_bonds(*mol_, unk_indices);
      if (should_initialize_ringinfo(new_mol.getNumAtoms())) {
        auto vec = match_atom_types(new_mol);
        for (size_t i = 0; i < unk_indices.size(); ++i) {
          atom_types.add_data(*mol_, mol_->getAtomWithIdx(unk_indices[i]), vec[i]);
        }
      }
    }

    rings_vec = populate_ring_entities();
  }

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
