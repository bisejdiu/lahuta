#ifndef LAHUTA_HPP
#define LAHUTA_HPP

#include <memory>
#include <string>
#include <vector>

#include <gemmi/mmread_gz.hpp>
#include <gemmi/model.hpp>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/BondIterators.h>

#include "GraphMol/Rings.h"
#include "aromatics.hpp"
#include "bond_order.hpp"
#include "bonds.hpp"
#include "contacts/features.hpp"
#include "contacts/interactions.hpp"
#include "convert.hpp"
#include "neighbors.hpp"
#include "nsgrid.hpp"
#include "ob/clean_mol.hpp"
#include "ob/kekulize.h"
#include "rings.hpp"

// selection parser
#include "visitor.hpp"

// contacts
#include "contacts/contacts.hpp"
#include "contacts/groups.hpp"

#include "nn.hpp"

#define LAHUTA_VERSION "0.14.0"
#define t() std::chrono::high_resolution_clock::now()
#define to_ms(d) std::chrono::duration_cast<std::chrono::milliseconds>(d)

namespace lahuta {

static float BONDED_NS_CUTOFF = 4.5;

class Topology {
public:
  std::vector<AtomType> atom_types;
  RingDataVec rings_vec;
  const RDKit::RWMol *mol;
  const Residues *residues;

  Topology(const RDKit::RWMol &mol) : mol(&mol) {}
  Topology() = default;
  ~Topology() { delete residues;}

public:
  void build_residues(const RDKit::RWMol &mol) {
    residues = new Residues(mol);
  }

  void assign_arpeggio_atom_types(RDKit::RWMol &mol) {

    std::vector<AtomType> indices(mol.getNumAtoms(), AtomType::NONE);
    std::vector<int> invalid_indices;
    Rings rings;

    // First pass: populate `indices` and track invalid & non-hydrogen atoms.
    for (auto atom : mol.atoms()) {
      // NOTE: `get_atom_type` expects standard residues
      AtomType atom_type = get_atom_type(atom);
      indices[atom->getIdx()] = atom_type;
      if (atom_type == AtomType::INVALID && atom->getAtomicNum() != 1) {
        invalid_indices.push_back(atom->getIdx());
      }

      if (AtomTypeFlags::has(atom_type, AtomType::AROMATIC)) {
        rings.add_ring_atom(atom);
      }
    }
    rings.process_rings(mol);
    RingDataVec _rings_vec = rings.get_rings_vector();

    // process invalid indices using a new RDKit object.
    if (!invalid_indices.empty()) {
      auto new_mol = filter_with_bonds(mol, invalid_indices);
      if (!new_mol.getRingInfo()->isInitialized()) {
        RDKit::MolOps::symmetrizeSSSR(new_mol);
      }

      auto vec = match_atom_types(new_mol);

      for (size_t i = 0; i < invalid_indices.size(); ++i) {
        indices[invalid_indices[i]] = vec[i];
      }

      auto new_mol_rings = new_mol.getRingInfo()->atomRings();
      for (auto &ring : new_mol_rings) {
        RDGeom::Point3D center, norm;
        Rings::compute_center(&new_mol, ring, center);
        Rings::compute_normal(&new_mol, ring, center, norm);
        std::vector<int> mapped_ring;
        std::vector<const RDKit::Atom *> atoms;
        for (const int &atom_idx : ring) {
          auto mappped_idx = invalid_indices[atom_idx];
          mapped_ring.push_back(mappped_idx);
          atoms.push_back(mol.getAtomWithIdx(mappped_idx));
        }

        /*RingData ring_data{center, norm, mapped_ring};*/
        RingData ring_data{center, norm, atoms};
        _rings_vec.rings.push_back(ring_data);
      }
    }

    atom_types = std::move(indices);
    rings_vec = std::move(_rings_vec);
  }

  void assign_molstar_typing() {

    // FIX: to be replaced by the entitytype manager
    std::vector<AtomType> new_atom_types = AtomTypeAnalysis::analyze(*mol);
    FeatureVec features = GroupTypeAnalysis::analyze(*mol, *residues);

    // populate rings_vec
    RingDataVec rings;
    for (const auto &ring : mol->getRingInfo()->atomRings()) {
      RDGeom::Point3D center, norm;
      Rings::compute_center(mol, ring, center);
      Rings::compute_normal(mol, ring, center, norm);

      std::vector<const RDKit::Atom *> atoms;
      for (const int &atom_idx : ring) {
        atoms.push_back(mol->getAtomWithIdx(atom_idx));
      }

      rings.rings.emplace_back(center, norm, atoms);
    }

    atom_types = std::move(new_atom_types);
    rings_vec = std::move(rings);
  }

  static void compute_bonds(RDKit::RWMol &mol, const NSResults &neighborResults) {
    auto result = assign_bonds(mol, neighborResults);

    // FIX: Refactor!

    // for (gemmi::Connection &conn : st.connections) {
    //   // Iterate over all models
    //   gemmi::Atom *a1 = st.first_model().find_cra(conn.partner1).atom;
    //   gemmi::Atom *a2 = st.first_model().find_cra(conn.partner2).atom;
    //
    //   if (mol.getBondBetweenAtoms(a1->serial-1, a2->serial-1) == nullptr) {
    //     mol.addBond((unsigned int)a1->serial-1, (unsigned int)a2->serial-1,
    //                 RDKit::Bond::BondType::SINGLE);
    //   }
    // }

    mol.updatePropertyCache(false);
    /*clean_bonds(mol, mol.getConformer());*/
    MolOps::setHybridization(mol);
    cleanup_predef(mol);

    if (result.has_unlisted_resnames) {
      clean_bonds(result.mol, result.mol.getConformer());
      perceive_bond_orders_obabel(result.mol);
      cleanup(result.mol);

      /*if (!result.mol.getRingInfo()->isInitialized()) {*/
      /*  std::cout << "Symmetrize SSSR\n";*/
      /*  RDKit::MolOps::symmetrizeSSSR(result.mol);*/
      /*}*/
      /*auto new_mol_rings = result.mol.getRingInfo()->atomRings();*/
      /*std::cout << "OB+RDKit HEME : \n";*/
      /*for (const auto &ring : new_mol_rings) {*/
      /*  std::cout << "";*/
      /*  for (const int &atom_idx : ring) {*/
      /*    std::cout << atom_idx << " ";*/
      /*  }*/
      /*  std::cout << "\n";*/
      /*  std::cout << " ";*/
      /*  for (const int &atom_idx : ring) {*/
      /*    std::cout << result.mol.getAtomWithIdx(atom_idx)->getIsAromatic() << " ";*/
      /*  }*/
      /*  std::cout << "\n";*/
      /*}*/

      merge_bonds(mol, result.mol, result.atom_indices);
    }
  }

private:
  static void cleanup_predef(RWMol &mol) {
    ROMol::VERTEX_ITER atBegin, atEnd;
    boost::tie(atBegin, atEnd) = mol.getVertices();
    while (atBegin != atEnd) {
      RDKit::Atom *atom = mol[*atBegin];
      atom->calcExplicitValence(false);

      // correct four-valent neutral N -> N+
      // This was github #1029
      if (atom->getAtomicNum() == 7 && atom->getFormalCharge() == 0 && atom->getExplicitValence() == 4) {
        atom->setFormalCharge(1);
      }
      ++atBegin;
    }
  }

  static void cleanup(RDKit::RWMol &mol) {
    cleanup_predef(mol);

    // FIXME: provide the mechanism as an option?
    // MolOps::fastFindRings(mol);
    // MolOps::findSSSR(mol);
    bool include_dative_bonds = true;
    MolOps::symmetrizeSSSR(mol, include_dative_bonds);
    MolOps::setAromaticity(mol);
  }

  static void
  merge_bonds(RDKit::RWMol &targetMol, RDKit::RWMol &sourceMol, const std::vector<int> &indexMap) {
    for (auto bondIt = sourceMol.beginBonds(); bondIt != sourceMol.endBonds(); ++bondIt) {
      const RDKit::Bond *bond = *bondIt;
      int bIdx = indexMap[bond->getBeginAtomIdx()];
      int eIdx = indexMap[bond->getEndAtomIdx()];
      if (targetMol.getBondBetweenAtoms(bIdx, eIdx) == nullptr) {

        // FIX: todo: instead of doing this check, we should iterate once over
        // all atoms and set the number of explicit hydrogens based on the
        // number of explicitly bonded hydrogens
        auto a = targetMol.getAtomWithIdx(bIdx);
        auto b = targetMol.getAtomWithIdx(eIdx);
        int is_a_h = a->getAtomicNum() == 1;
        int is_b_h = b->getAtomicNum() == 1;

        if (is_a_h ^ is_b_h) {
          auto non_h_atom = a->getAtomicNum() == 1 ? b : a;
          non_h_atom->setNumExplicitHs(non_h_atom->getNumExplicitHs() + 1);
        }
        targetMol.addBond(bIdx, eIdx, bond->getBondType());
      }
    }
  }
};

class Luni {
  // FIX: move down
private:
  std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();
  Structure st;
  NSResults bonded_nps;
  double _cutoff;
  FastNS grid; // FIXME: is this needed?
  Topology topology = Topology(*mol);
  /*std::vector<Feature> features;*/
  FeatureVec features;

  void process_file(std::string file_path) {
    RDKit::Conformer *conformer = new RDKit::Conformer();
    auto start = std::chrono::high_resolution_clock::now();
    st = read_structure_gz(file_path);
    file_name = file_path;
    std::cout << "Read Structure using gemmi: " << to_ms(t() - start).count() << "\n";

    start = std::chrono::high_resolution_clock::now();
    gemmiStructureToRDKit(*mol, st, *conformer, false);
    std::cout << "Convert gemmi to RDKit: " << to_ms(t() - start).count() << "\n";

    start = std::chrono::high_resolution_clock::now();
    mol->updatePropertyCache(false);
    std::cout << "Update Property Cache: " << to_ms(t() - start).count() << "\n";

    start = std::chrono::high_resolution_clock::now();
    mol->addConformer(conformer, true);
    std::cout << "Add Conformer: " << to_ms(t() - start).count() << "\n";
  }

  void create_topology() {

    auto start = std::chrono::high_resolution_clock::now();
    const auto &conf = get_conformer();
    grid = FastNS(conf.getPositions(), _cutoff);
    bonded_nps = grid.self_search();
    std::cout << "Self Search: " << to_ms(t() - start).count() << "\n";

    start = std::chrono::high_resolution_clock::now();
    Topology::compute_bonds(*mol, bonded_nps);
    std::cout << "Compute Bonds: " << to_ms(t() - start).count() << "\n";



    topology.build_residues(*mol);
    initialize_and_populate_ringinfo(*mol, *topology.residues);


    /*if (!mol->getRingInfo()->isInitialized()) {*/
    /*  RDKit::MolOps::symmetrizeSSSR(*mol);*/
    /*}*/
    // print all rings
    /*std::cout << "x--> Rings: " << mol->getRingInfo()->numRings() << "\n";*/
    /*for (const auto &ring : mol->getRingInfo()->atomRings()) {*/
    /*  std::cout << "x--> Ring: ";*/
    /*  auto first = ring.front();*/
    /*  auto first_atom = mol->getAtomWithIdx(first);*/
    /*  auto info = static_cast<RDKit::AtomPDBResidueInfo *>(first_atom->getMonomerInfo());*/
    /*  std::cout << info->getResidueName() << " " << info->getResidueNumber() << " ";*/
    /**/
    /*  for (const auto &atom_idx : ring) {*/
    /*    auto atom = mol->getAtomWithIdx(atom_idx);*/
    /*    std::cout << atom->getIdx() << " ";*/
    /*  }*/
    /*  std::cout << "\n";*/
    /*}*/

    start = std::chrono::high_resolution_clock::now();
    /*topology.assign_atom_types(*mol);*/

    topology.assign_molstar_typing();
    /*topology.assign_arpeggio_atom_types(*mol);*/

    std::cout << "Assign Atom Types: " << to_ms(t() - start).count() << "\n";

    for (const auto &atom : mol->atoms()) {
      /*auto is_conj = isConjugated(*mol, *atom);*/
      /*std::cout << "Conjugated: " << is_conj << "\n";*/
      auto info = static_cast<RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());

      /*auto no_ih = calculateHydrogensCharge(*mol, *atom, true, true);*/
      /*calculateHydrogensCharge(*mol, *atom, true, true);*/
      /*auto h_count = get_h_count(*mol, *atom);*/
      /*std::cout << info->getName() << " " << */
      /*  h_count << " ? " <<*/
      /*  no_ih << " ---- " <<*/
      /*  no_ih + atom->getNumExplicitHs() << " | " <<*/
      /*  atom->getNumImplicitHs() << " " <<*/
      /*  atom->getNumExplicitHs() << " " <<*/
      /*  atom->getTotalNumHs() << " " << */
      /*  "\n";*/
    }
  }

public:
  Luni() = default; // FIX: remove?
  explicit Luni(std::string file_name) : _cutoff(BONDED_NS_CUTOFF) {
    process_file(file_name);

    auto start = std::chrono::high_resolution_clock::now();
    create_topology();

    // FIX: double call to Residues(*mol)
    std::cout << "Create Topology: " << to_ms(t() - start).count() << "\n";
    Residues residues(*mol);
    features = std::move(GroupTypeAnalysis::analyze(*mol, residues));
    std::cout << "Analyze Group Types: " << to_ms(t() - start).count() << "\n";
    std::cout << "Features: " << get_features().features.size() << "\n";
    /*features = GroupTypeAnalyzer::*/
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
  const RingDataVec &get_rings() const { return topology.rings_vec; }

  /*const RingDataVec new_get_rings() const {*/
  /*  RingDataVec rings;*/
  /*  for (const auto &ring : mol->getRingInfo()->atomRings()) {*/
  /*    RDGeom::Point3D center, norm;*/
  /*    Rings::compute_center(mol.get(), ring, center);*/
  /*    Rings::compute_normal(mol.get(), ring, center, norm);*/
  /**/
  /*    std::vector<const RDKit::Atom *> atoms;*/
  /*    for (const int &atom_idx : ring) {*/
  /*      atoms.push_back(mol->getAtomWithIdx(atom_idx));*/
  /*    }*/
  /**/
  /*    rings.rings.emplace_back(center, norm, atoms);*/
  /*  }*/
  /*  return rings;*/
  /*};*/

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
    auto centers = rings.centers();

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
    auto centers = rings.centers();

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
  const FeatureVec &get_features() const { return features; }

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

    FeatureVec group_features = GroupTypeAnalysis::analyze(*mol, Residues(*mol));
    Interactions processor(this, InteractionOptions{5.0});
    auto ionic = processor.find_ionic_interactions();
    auto hbonds = processor.find_hbond_interactions();
    c.add(ionic);
    c.add(hbonds);

    return c;
  }

  /*void assign_molstar_atom_types() { topology.assign_molstar_atom_types(*mol); }*/
  void assign_molstar_atom_types() { topology.assign_molstar_typing(); }
  void assign_arpeggio_atom_types() { topology.assign_arpeggio_atom_types(*mol); }

private:
  std::vector<int> filtered_indices;
  std::unordered_map<lahuta::EntityType, std::vector<EntityID>> entities;
};

} // namespace lahuta

#endif // LAHUTA_HPP
