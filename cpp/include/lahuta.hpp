#ifndef LAHUTA_HPP
#define LAHUTA_HPP

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <gemmi/mmread_gz.hpp>
#include <gemmi/model.hpp>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/BondIterators.h>

#include "convert.hpp"
#include "neighbors.hpp"
#include "topology.hpp"
#include "spdlog/spdlog.h"

#define LAHUTA_VERSION "0.15.0"

namespace lahuta {

class Luni { // rename to Lahuta
public:
  explicit Luni(std::string file_name) : file_name_(file_name) {
    spdlog::info("Processing file: {}", file_name_);
    read_structure();
  }

  bool build_topology(std::optional<TopologyBuildingOptions> tops = TopologyBuildingOptions()) {
    try {
      topology.emplace(mol);
      topology->build(tops.value_or(TopologyBuildingOptions()));
    } catch (const std::runtime_error &e) {
      return false;
    }
    return true;
  }

  static Luni create(std::shared_ptr<RDKit::RWMol> mol) {
    return Luni(mol);
  }

  static Luni create(const gemmi::Structure &st) {
    auto mol = std::make_shared<RDKit::RWMol>();
    RDKit::Conformer *conformer = new RDKit::Conformer();
    create_RDKit_repr(*mol, st, *conformer, false);
    mol->updatePropertyCache(false);
    mol->addConformer(conformer, true);
    return Luni(mol);
  }

  static Luni create(const IR &ir) {
    auto mol = std::make_shared<RDKit::RWMol>();
    IR_to_RWMol(*mol, ir);
    return Luni(mol);
  }

  std::string get_file_name() { return file_name_; };
  const Topology &get_topology() const { return *get_topology_ptr(); }
  bool has_topology_built() { return topology.has_value(); };

  RDKit::RWMol       &get_molecule()       { return *mol; }
  const RDKit::RWMol &get_molecule() const { return *mol; }

  RDKit::Conformer       &get_conformer(int id = -1)       { return mol->getConformer(id); }
  const RDKit::Conformer &get_conformer(int id = -1) const { return mol->getConformer(id); }

  const AtomEntityCollection  &get_atom_types() const { return get_topology_ptr()->get_atom_types(); }
  const RingEntityCollection  &get_rings()      const { return get_topology_ptr()->get_rings(); }
  const GroupEntityCollection &get_features()   const { return get_topology_ptr()->get_features(); }

  // FIX: should probably get rid of these
  Neighbors<AtomAtomPair> find_neighbors(double cutoff, int res_dif) {
    NSResults ns = find_neighbors_opt(cutoff);
    if (res_dif > 0) {
      ns = remove_adjascent_residueid_pairs(ns, res_dif);
    }
    return Neighbors<AtomAtomPair>(*this, std::move(ns), false);
  }

  NSResults find_neighbors2(double cutoff, int res_dif) {
    NSResults ns = find_neighbors_opt(cutoff);
    if (res_dif > 0) {
      ns = remove_adjascent_residueid_pairs(ns, res_dif);
    }
    return ns;
  }

  // FIX: move this to a separate class or namespace

  //! Returns the atoms of the molecule.
  const auto atoms() const { return mol->atoms(); }

  //! Returns the number of atoms in the molecule.
  const auto n_atoms() const { return mol->getNumAtoms(); }

  //! Returns the names of the atoms.
  const std::vector<std::string> names() const;

  //! Returns the symbols of the atoms.
  const std::vector<std::string> symbols() const;

  //! Returns the residue indices of the atoms.
  const std::vector<int> indices() const;

  //! Returns the atomic numbers of the atoms.
  const std::vector<int> atomic_numbers() const;

  //! Returns the elements of the atoms.
  const std::vector<std::string> elements() const;

  //! Returns the residue names of the atoms.
  const std::vector<std::string> resnames() const;

  //! Returns the residue ids of the atoms.
  const std::vector<int> resids() const;

  //! Returns the residue indices of the atoms.
  const std::vector<int> resindices() const;

  //! Returns the chain labels of the atoms.
  const std::vector<std::string> chainlabels() const;

  const std::vector<RDGeom::Point3D> &positions(int confId = -1) const {
    return get_conformer(confId).getPositions();
  }

  //! get the total number of bonded h atoms (explicit + implicit) for all atoms in the system
  std::vector<int> total_hydrogen_count() const {
    std::vector<int> hydrogen_counts;
    for (const auto &atom : mol->atoms()) {
      hydrogen_counts.push_back(atom->getTotalNumHs());
    }
    return hydrogen_counts;
  }

  template <typename T> friend class Neighbors;
  friend class Contacts;

  // FIX: add helper functions to get topology information
  // FIX: Move these to the topology class
  const RDKit::Atom *get_atom(int idx) const { return mol->getAtomWithIdx(idx); }

private:
  explicit Luni(std::shared_ptr<RDKit::RWMol> valid_mol) : mol(valid_mol), topology(valid_mol) {}

  template <typename T> std::vector<T> atom_attrs(std::function<T(const RDKit::Atom *)> func) const;

  template <typename T>
  std::vector<std::reference_wrapper<const T>>
  atom_attrs_ref(std::function<const T &(const RDKit::Atom *)> func) const;

  // FIX: make find_neighbors_opt public?
  NSResults find_neighbors_opt(double cutoff = BONDED_NEIGHBOR_SEARCH_CUTOFF);
  NSResults remove_adjascent_residueid_pairs(NSResults &results, int res_diff);

  auto match_smarts_string(std::string sm, std::string atype = "", bool log_values = false) const;

  const Topology* get_topology_ptr() const {
    if (!topology) {
      spdlog::error("Cannot return topology, because no topology has been built.");
      throw std::runtime_error("Topology not built");
    }
    return &topology.value();
  }

  void read_structure() {
    auto st = gemmi::read_structure_gz(file_name_);

    RDKit::Conformer *conformer = new RDKit::Conformer();
    create_RDKit_repr(*mol, st, *conformer, false);
    mol->updatePropertyCache(false);
    mol->addConformer(conformer, true);
  }

public:
  std::string file_name_;

  // FIX: document these four
  static std::vector<std::string> tokenize(const std::string &str);
  static std::vector<std::string> tokenize_simple(const std::string &str);
  std::vector<int> parse_and_filter(const std::string &selection) const;
  bool parse_expression(const std::string &selection);

  // FIX: document the difference
  Luni filter() const;
  Luni filter_luni(const std::vector<int> &atom_indices) const;

  static std::vector<int> factorize(const std::vector<std::string> &labels);
  static std::vector<std::string> find_elements(const std::vector<int> &atomic_numbers);

  static int count_unique(const std::vector<int> &vec);
  static int count_unique(const std::vector<std::string> &vec);

  // FIX: what do we do about the entities here?
  template <typename T> const T &get_entity(EntityID id) const;
  const std::vector<EntityID> &get_atom_entities();
  const std::vector<EntityID> &get_ring_entities();
  const std::vector<EntityID> &get_group_entities();

  // FIX: part of the topology class. Debatable if they should be here. Can simply be called via the topology attribute
  void assign_molstar_atom_types()  { 
    if (topology) { topology->assign_molstar_typing(); } 
    else { spdlog::error("Topology not built. Cannot assign Molstar atom types."); }
  }

  void assign_arpeggio_atom_types() {
    if (topology) { topology->assign_arpeggio_atom_types(); } 
    else { spdlog::error("Topology not built. Cannot assign Arpeggio atom types."); }
  }

private:
  std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();
  std::optional<Topology> topology;
  std::vector<int> filtered_indices;
  std::unordered_map<lahuta::EntityType, std::vector<EntityID>> entities;
};

} // namespace lahuta

#endif // LAHUTA_HPP
