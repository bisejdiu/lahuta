#ifndef LAHUTA_HPP
#define LAHUTA_HPP

#include <memory>
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

class Luni {
  // FIX: move down
private:
  std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();
  std::shared_ptr<NSResults> neighbors;

  double _cutoff;
  FastNS grid; // FIXME: is this needed?
  Topology topology;
  GroupEntityCollection features;

  void process_file(std::string file_path_) {
    file_name = file_path_;
    auto st = gemmi::read_structure_gz(file_path_);

    RDKit::Conformer *conformer = new RDKit::Conformer();
    create_RDKit_repr(*mol, st, *conformer, false);
    mol->updatePropertyCache(false);
    mol->addConformer(conformer, true);

    topology = Topology(mol.get());
  }

  void create_topology() {

    try {
      grid = FastNS(get_conformer().getPositions(), _cutoff);
      neighbors = std::make_shared<NSResults>(grid.self_search());

      Topology::compute_bonds(*mol, *neighbors);

      topology.build_residues(*mol);
      initialize_and_populate_ringinfo(*mol, *topology.residues);

      topology.assign_molstar_typing();
      /*topology.assign_arpeggio_atom_types();*/

    } catch (const std::runtime_error &e) {
      throw std::runtime_error("Error creating topology: " + std::string(e.what()));
    }
  }

public:
  Luni() = default; // FIX: remove?
  explicit Luni(std::string file_name) : _cutoff(BONDED_NS_CUTOFF) {
    spdlog::info("Processing file: {}", file_name);
    process_file(file_name);

    try {
      create_topology();
    } catch (const std::runtime_error &e) {
      std::cerr << "Critical error: " << e.what() << "\n";
      success = false;
      return;
    }

    // FIX: double call to Residues(*mol)
    Residues residues(*mol);
    features = std::move(GroupTypeAnalysis::analyze(*mol, residues));
  }

  // Luni from IR:
  Luni(const IR &ir) : _cutoff(BONDED_NS_CUTOFF) {
    IR_to_RWMol(*mol, ir);
    create_topology();
  }

  //! get the total number of bonded h atoms (explicit + implicit) for all atoms
  //! in the system
  std::vector<int> total_hydrogen_count() const {
    std::vector<int> hydrogen_counts;
    for (const auto &atom : mol->atoms()) {
      hydrogen_counts.push_back(atom->getTotalNumHs());
    }
    return hydrogen_counts;
  }

  const std::vector<RDGeom::Point3D> &positions(int confId = -1) const {
    return get_conformer(confId).getPositions();
  }

  const RDKit::Conformer &get_conformer(int id = -1) const { return mol->getConformer(id); }

  RDKit::RWMol &get_molecule() { return *mol; }
  const RDKit::RWMol &get_molecule() const { return *mol; }

  const double get_cutoff() const { return _cutoff; }

  const std::vector<AtomType> &get_atom_types() const { return topology.atom_types; }
  const RingEntityCollection &get_rings() const { return topology.rings_vec; }

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

  // FIX: computed distances represent atom-ring center distances.
  Neighbors<AtomRingPair> find_ring_neighbors(double cutoff, int res_dif = 1) {

    auto rings = topology.rings_vec;
    // FIX: keep an instance of the grid in the class
    auto grid = FastNS(mol->getConformer().getPositions(), cutoff);
    auto centers = rings.positions();

    NSResults nbrs = grid.search(centers);
    if (res_dif > 0) {
      nbrs = remove_adjascent_residueid_pairs(nbrs, res_dif);
    }
    return Neighbors<AtomRingPair>(*this, std::move(nbrs), false);
  }

  // FIX: computed distances represent atom-ring center distances.
  // FIX: neighbors contain atoms that are also part of the ring.
  auto find_ring_neighbors2(double cutoff, int res_dif = 1) {

    auto rings = topology.rings_vec;
    // FIX: keep an instance of the grid in the class
    auto grid = FastNS(mol->getConformer().getPositions(), cutoff);
    auto centers = rings.positions();

    NSResults nbrs = grid.search(centers);
    // log the first 10 pairs of neighbors if that many exist
    for (size_t i = 0; i < nbrs.size(); ++i) {
      auto &pair = nbrs.get_pairs()[i];
      /*std::cout << "Ring: " << pair.first << " Atom: " << pair.second << "\n";*/
      /*auto atom = mol->getAtomWithIdx(pair.first);*/
      /*auto &ring = rings.rings[pair.second];*/
      /*std::cout << "Atom: " << atom->getIdx() << " Ring: " << ring.center << "\n";*/
    }
    if (res_dif > 0) {
      nbrs = remove_adjascent_residueid_pairs(nbrs, res_dif);
    }
    return nbrs;
  }

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

  template <typename T> friend class Neighbors;
  friend class Contacts;
  /*friend class InteractionContainer::test_interface;*/

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
  std::string file_name;
  static std::vector<std::string> tokenize(const std::string &str);
  static std::vector<std::string> tokenize_simple(const std::string &str);
  std::vector<int> parse_and_filter(const std::string &selection) const;
  bool parse_expression(const std::string &selection);
  Luni filter_luni(const std::vector<int> &atom_indices) const;
  Luni filter() const;
  static std::vector<int> factorize(const std::vector<std::string> &labels);
  static int count_unique(const std::vector<int> &vec);
  static int count_unique(const std::vector<std::string> &vec);
  static std::vector<std::string> find_elements(const std::vector<int> &atomic_numbers);

  /*const std::vector<Feature> &get_features() const { return features; }*/
  const GroupEntityCollection &get_features() const { return features; }

  template <typename T> const T &get_entity(EntityID id) const;
  const std::vector<EntityID> &get_atom_entities();
  const std::vector<EntityID> &get_ring_entities();
  const std::vector<EntityID> &get_group_entities();

  Contacts test_find_neighbors(double cutoff, int res_dif) {
    Contacts c(this);
    /*const auto &atom_entities = get_atom_entities();*/
    /*auto atom_neighbors = find_neighbors2(6.0, 10);*/
    /*c.add_many(atom_neighbors, atom_entities);*/
    /**/
    /*std::vector<RingData> rings = get_rings().rings;*/
    /*const auto &ring_entities = get_ring_entities();*/
    /*auto ring_neighbors = find_ring_neighbors2(6.0);*/
    /*c.add_many(ring_neighbors, ring_entities, atom_entities);*/

    GroupEntityCollection group_features = GroupTypeAnalysis::analyze(*mol, Residues(*mol));
    Interactions processor(*this, InteractionOptions{5.0});
    auto ionic = processor.ionic();
    auto hbonds = processor.hbond();
    c.add(ionic);
    c.add(hbonds);

    return c;
  }

  /*void assign_molstar_atom_types() { topology.assign_molstar_atom_types(*mol); }*/
  void assign_molstar_atom_types() { topology.assign_molstar_typing(); }
  void assign_arpeggio_atom_types() { topology.assign_arpeggio_atom_types(); }

private:
  std::vector<int> filtered_indices;
  std::unordered_map<lahuta::EntityType, std::vector<EntityID>> entities;
};

} // namespace lahuta

#endif // LAHUTA_HPP
