#ifndef LAHUTA_HPP
#define LAHUTA_HPP

#include <memory>
#include <string>
#include <vector>

#include <gemmi/mmread_gz.hpp> // for read_structure_gz
#include <gemmi/model.hpp>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/BondIterators.h>

#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/join.hpp>

#include "bond_order.hpp"
#include "bonds.hpp"
#include "convert.hpp"
#include "neighbors.hpp"
#include "nsgrid.hpp"
#include "ob/clean_mol.hpp"
#include "ob/kekulize.h"
#include "rings.hpp"

// selection parser
#include "parser.hpp"
#include "visitor.hpp"

#define LAHUTA_VERSION "0.11.0"
#define t() std::chrono::high_resolution_clock::now()
#define to_ms(d) std::chrono::duration_cast<std::chrono::milliseconds>(d)

namespace lahuta {

static float BONDED_NS_CUTOFF = 4.5;

class Topology {
public:
  std::vector<AtomType> atom_types;
  RingDataVec rings_vec;

  Topology() = default;
  Topology(const Topology &other) = default;
  Topology(Topology &&other) = default;
  Topology &operator=(const Topology &other) = default;
  Topology &operator=(Topology &&other) = default;
  ~Topology() = default;

public:
  void assign_atom_types(RDKit::RWMol &mol) {

    std::vector<AtomType> indices(mol.getNumAtoms(), AtomType::NONE);

    std::vector<int> invalid_indices;
    Rings rings;

    // First pass: populate `indices` and track invalid non-hydrogen atoms.
    for (auto atom : mol.atoms()) {
      AtomType atom_type = get_atom_type(atom);
      indices[atom->getIdx()] = atom_type;
      if (atom_type == AtomType::INVALID && atom->getAtomicNum() != 1) {
        invalid_indices.push_back(atom->getIdx());
      }

      if (AtomTypeFlags::has(atom_type, AtomType::AROMATIC)) {
        rings.add_ring_atom(atom);
        RDGeom::Point3D atom_pos =
            mol.getConformer().getAtomPos(atom->getIdx());
      }
    }
    rings.process_rings(mol);
    RingDataVec _rings_vec = rings.get_rings_vector();

    // Handle invalid indices using a new RDKit object.
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
        RDGeom::Point3D center, norm1, norm2;
        Rings::find_center_and_normal(new_mol.getConformer(), ring, center,
                                      norm1, norm2);
        std::vector<int> mapped_ring;
        for (const int &atom_idx : ring) {
          auto mappped_idx = invalid_indices[atom_idx];
          mapped_ring.push_back(mappped_idx);
        }

        RingData ring_data{center, norm1, norm2, mapped_ring};
        _rings_vec.rings.push_back(ring_data);
      }
    }

    // for (size_t i = 0; i < indices.size(); ++i) {
    // // for (size_t i = 0; i < 10; ++i) {
    //   auto atom = mol.getAtomWithIdx(i);
    //   auto *info =
    //       static_cast<RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    //
    //   if (atom->getAtomicNum() != 1) {
    //     std::cout << info->getResidueName() << " " << info->getName() << " "
    //               << atom_type_to_string(indices[i]) << std::endl;
    //   }
    // }

    atom_types = std::move(indices);
    rings_vec = std::move(_rings_vec);
  }

  static void compute_bonds(RDKit::RWMol &mol,
                            const NSResults &neighborResults) {
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
    cleanup_predef(mol);

    if (result.has_unlisted_resnames) {
      clean_bonds(result.mol, result.mol.getConformer());
      perceive_bond_orders_obabel(result.mol);
      cleanup(result.mol);

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
      if (atom->getAtomicNum() == 7 && atom->getFormalCharge() == 0 &&
          atom->getExplicitValence() == 4) {
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

  static void merge_bonds(RDKit::RWMol &targetMol, RDKit::RWMol &sourceMol,
                          const std::vector<int> &indexMap) {
    for (auto bondIt = sourceMol.beginBonds(); bondIt != sourceMol.endBonds();
         ++bondIt) {
      const RDKit::Bond *bond = *bondIt;
      int bIdx = indexMap[bond->getBeginAtomIdx()];
      int eIdx = indexMap[bond->getEndAtomIdx()];
      if (targetMol.getBondBetweenAtoms(bIdx, eIdx) == nullptr) {
        targetMol.addBond(bIdx, eIdx, bond->getBondType());
      }
    }
  }
};

class Luni {
private:
  std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();
  Structure st;
  NSResults bonded_nps;
  double _cutoff;
  FastNS grid; // FIXME: is this needed?
  Topology topology;

  void process_file(std::string file_name) {
    RDKit::Conformer *conformer = new RDKit::Conformer();
    auto start = std::chrono::high_resolution_clock::now();
    st = read_structure_gz(file_name);
    std::cout << "Read Structure using gemmi: " << to_ms(t() - start).count()
              << "\n";

    start = std::chrono::high_resolution_clock::now();
    gemmiStructureToRDKit(*mol, st, *conformer, false);
    std::cout << "Convert gemmi to RDKit: " << to_ms(t() - start).count()
              << "\n";

    start = std::chrono::high_resolution_clock::now();
    mol->updatePropertyCache(false);
    std::cout << "Update Property Cache: " << to_ms(t() - start).count()
              << "\n";

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

    start = std::chrono::high_resolution_clock::now();
    topology.assign_atom_types(*mol);
    std::cout << "Assign Atom Types: " << to_ms(t() - start).count() << "\n";
  }

public:
  Luni() = default;
  explicit Luni(std::string file_name) : _cutoff(BONDED_NS_CUTOFF) {
    process_file(file_name);

    auto start = std::chrono::high_resolution_clock::now();
    create_topology();
    std::cout << "Create Topology: " << to_ms(t() - start).count() << "\n";
  }

  const std::vector<RDGeom::Point3D> &positions(int confId = -1) const {
    return get_conformer(confId).getPositions();
  }

  const RDKit::Conformer &get_conformer(int id = -1) const {
    return mol->getConformer(id);
  }

  RDKit::RWMol &get_molecule() { return *mol; }
  const RDKit::RWMol &get_molecule() const { return *mol; }

  const double get_cutoff() const { return _cutoff; }

  const std::vector<AtomType> &get_atom_types() const {
    return topology.atom_types;
  }
  const RingDataVec &get_rings() const { return topology.rings_vec; }

  template <typename T>
  Neighbors<T> find_neighbors(double cutoff, int res_dif) {
    NSResults ns = find_neighbors_opt(cutoff);
    ns = remove_adjascent_residueid_pairs(ns, res_dif);
    return Neighbors<T>(*this, std::move(ns), false);
  }

  // FIX: computed distances represent atom-ring center distances.
  Neighbors<AtomRingPair> find_ring_neighbors(double cutoff, int res_dif=1) {

    auto rings = topology.rings_vec; 
    auto grid = FastNS(mol->getConformer().getPositions(), 6.0);
    auto centers = rings.centers_rkdit();

    NSResults nbrs = grid.search(centers);
    if (res_dif > 0) {
      nbrs = remove_adjascent_residueid_pairs(nbrs, res_dif);
    }
    return Neighbors<AtomRingPair>(*this, std::move(nbrs), false);
  }

  // FIX: move to source
  const RDKit::Atom &get_atom(int index) const {
    return *mol->getAtomWithIdx(index);
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

private:
  template <typename T>
  std::vector<T> atom_attrs(std::function<T(const RDKit::Atom *)> func) const;

  template <typename T>
  std::vector<std::reference_wrapper<const T>>
  atom_attrs_ref(std::function<const T &(const RDKit::Atom *)> func) const;

  // FIX: make find_neighbors_opt public?
  NSResults find_neighbors_opt(double cutoff = BONDED_NS_CUTOFF);
  NSResults remove_adjascent_residueid_pairs(NSResults &results, int res_diff);

  std::vector<RDKit::MatchVectType>
  match_smarts_string(std::string sm, std::string atype = "",
                      bool log_values = false) const;

public:
  static std::vector<std::string> tokenize(const std::string &str);
  static std::vector<std::string> tokenize_simple(const std::string &str);
  std::vector<int> parse_and_filter(const std::string &selection) const;
  bool parse_expression(const std::string &selection);
  Luni filter_luni(const std::vector<int> &atom_indices) const;
  Luni filter() const;

private:
  std::vector<int> filtered_indices;
};

} // namespace lahuta

#endif // LAHUTA_HPP
