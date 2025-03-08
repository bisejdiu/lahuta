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

#include "aromatics.hpp"
#include "contacts/interactions.hpp"
#include "convert.hpp"
#include "neighbors.hpp"
#include "nsgrid.hpp"
#include "ob/kekulize.h"

#include "contacts/groups.hpp"
#include "nn.hpp"
#include "spdlog/spdlog.h"
#include "topology.hpp"
#include "visitor.hpp" // selection parser (bad file name)

#define LAHUTA_VERSION "0.15.0"

namespace lahuta {

static float BONDED_NS_CUTOFF = 4.5;

enum ContactComputerType { Arpeggio, Molstar };

class Luni {
  // FIX: move down
private:
  void read_structure() {
    auto st = gemmi::read_structure_gz(file_name_);

    RDKit::Conformer *conformer = new RDKit::Conformer();
    create_RDKit_repr(*mol, st, *conformer, false);
    mol->updatePropertyCache(false);
    mol->addConformer(conformer, true);

    topology = Topology(mol.get()); // initialize the topology (not built yet)
  }

  void create_topology() {

    try {
      // NOTE: GRID
      grid = FastNS(get_conformer().getPositions());
      auto ok = grid.build(_cutoff);
      if (!ok) {
          throw std::runtime_error("Box dimension too small for the given cutoff.");
      }
      // FIX: perhaps we should not store neighbors (large size, memory usage, etc.)
      // Further, neighbors are only used in bond computation
      neighbors = std::make_shared<NSResults>(grid.self_search());

      // NOTE: BONDS
      Topology::compute_bonds(*mol, *neighbors);

      // NOTE: RESIDUES
      topology.build_residues(*mol);

      // NOTE: RINGS
      initialize_and_populate_ringinfo(*mol, *topology.residues);

      // NOTE: ATOM TYPES
      if (c_type == ContactComputerType::Arpeggio) {
        topology.assign_arpeggio_atom_types();
      } else {
        topology.assign_molstar_typing();
      }

    } catch (const std::runtime_error &e) {
      throw std::runtime_error("Error creating topology: " + std::string(e.what()));
    }
  }

public:
  Luni() : _cutoff(BONDED_NS_CUTOFF) {}
  // FIX: should to `build` the topology same as I am building the grid
  explicit Luni(std::string file_name, std::optional<ContactComputerType> c_type_ = std::nullopt)
    : _cutoff(BONDED_NS_CUTOFF), c_type(c_type_.value_or(ContactComputerType::Molstar)), file_name_(file_name) {

    spdlog::info("Processing file: {}", file_name_);
    read_structure();

    try {
      create_topology();
    } catch (const std::runtime_error &e) {
      std::cerr << "Critical error: " << e.what() << "\n";
      success = false;
      return;
    }

    // FIX: double call to Residues(*mol)
    Residues residues(*mol);
    if (!residues.build()) {
          throw std::runtime_error("Could not build residue information!");
    }
    features = std::move(GroupTypeAnalysis::analyze(*mol, residues));
  }

  static Luni build(std::shared_ptr<RDKit::RWMol> mol) {
    Luni luni;

    luni.mol = mol;
    luni.topology = Topology(luni.mol.get());

    luni.create_topology();
    return luni;
  }

  static Luni build(const gemmi::Structure &st) {
    Luni luni;

    RDKit::Conformer *conformer = new RDKit::Conformer();
    create_RDKit_repr(*luni.mol, st, *conformer, false);
    luni.mol->updatePropertyCache(false);
    luni.mol->addConformer(conformer, true);
    luni.topology = Topology(luni.mol.get());

    luni.create_topology();

    Residues residues(*luni.mol);
    if (!residues.build()) {
          throw std::runtime_error("Could not build residue information!");
    }
    luni.features = std::move(GroupTypeAnalysis::analyze(*luni.mol, residues));

    return luni;
  }

  // Luni from IR:
  Luni(const IR &ir) : _cutoff(BONDED_NS_CUTOFF) {
    IR_to_RWMol(*mol, ir);
    create_topology();
  }

  Topology       &get_topology()       { return topology; }
  const Topology &get_topology() const { return topology; };

  RDKit::RWMol       &get_molecule()       { return *mol; }
  const RDKit::RWMol &get_molecule() const { return *mol; }

  RDKit::Conformer       &get_conformer(int id = -1)       { return mol->getConformer(id); }
  const RDKit::Conformer &get_conformer(int id = -1) const { return mol->getConformer(id); }

  const double get_cutoff() const { return _cutoff; }

  const std::vector<AtomType> &get_atom_types() const { return topology.atom_types; }
  const RingEntityCollection  &get_rings()      const { return topology.rings_vec; }

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
  template <typename T> std::vector<T> atom_attrs(std::function<T(const RDKit::Atom *)> func) const;

  template <typename T>
  std::vector<std::reference_wrapper<const T>>
  atom_attrs_ref(std::function<const T &(const RDKit::Atom *)> func) const;

  // FIX: make find_neighbors_opt public?
  NSResults find_neighbors_opt(double cutoff = BONDED_NS_CUTOFF);
  NSResults remove_adjascent_residueid_pairs(NSResults &results, int res_diff);

  auto match_smarts_string(std::string sm, std::string atype = "", bool log_values = false) const;

public:
  bool success{true};
  std::string file_name_;

  static std::vector<std::string> tokenize(const std::string &str);
  static std::vector<std::string> tokenize_simple(const std::string &str);

  std::vector<int> parse_and_filter(const std::string &selection) const;
  bool parse_expression(const std::string &selection);

  Luni filter_luni(const std::vector<int> &atom_indices) const;
  Luni filter() const;

  static std::vector<int> factorize(const std::vector<std::string> &labels);
  static std::vector<std::string> find_elements(const std::vector<int> &atomic_numbers);

  static int count_unique(const std::vector<int> &vec);
  static int count_unique(const std::vector<std::string> &vec);

  const GroupEntityCollection &get_features() const { return features; }

  template <typename T> const T &get_entity(EntityID id) const;
  const std::vector<EntityID> &get_atom_entities();
  const std::vector<EntityID> &get_ring_entities();
  const std::vector<EntityID> &get_group_entities();

  void assign_molstar_atom_types()  { topology.assign_molstar_typing(); }
  void assign_arpeggio_atom_types() { topology.assign_arpeggio_atom_types(); }

private:
  std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();
  std::shared_ptr<NSResults> neighbors;

  double _cutoff;
  ContactComputerType c_type;
  FastNS grid; // FIXME: is this needed?
  Topology topology;
  GroupEntityCollection features;

  std::vector<int> filtered_indices;
  std::unordered_map<lahuta::EntityType, std::vector<EntityID>> entities;
};

} // namespace lahuta

#endif // LAHUTA_HPP
