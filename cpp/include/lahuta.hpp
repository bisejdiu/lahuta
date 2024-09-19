#ifndef LAHUTA_HPP
#define LAHUTA_HPP

#include <chrono>
#include <memory>
#include <string>
#include <vector>

#include <gemmi/mmread_gz.hpp> // for read_structure_gz
#include <rdkit/GraphMol/BondIterators.h>

#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/join.hpp>

#include "Geometry/point.h"
#include "bond_order.hpp"
#include "bonds.hpp"
#include "convert.hpp"
#include "gemmi/model.hpp"
#include "gemmi/seqid.hpp"
#include "nsgrid.hpp"
#include "rings.hpp"
#include "neighbors.hpp"
#include "ob/clean_mol.hpp"
#include "ob/kekulize.h"

#define LAHUTA_VERSION "0.11.0"

namespace lahuta {

static float BONDED_NS_CUTOFF = 4.5;

inline void basic_mol_cleanup(RWMol &mol) {
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

class Topology {
public:
  std::vector<AtomType> atom_types;
  RingDataVec rings_vec;

public:
  static void cleanup_predef(RDKit::RWMol &mol) { basic_mol_cleanup(mol); }

  static void cleanup(RDKit::RWMol &mol) {
    basic_mol_cleanup(mol);

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

      // if (has(atom_type, AtomType::CYCLICAL) || has(atom_type, AtomType::AROMATIC)) {
      if (AtomTypeFlags::has(atom_type, AtomType::AROMATIC)) {
        rings.add_ring_atom(atom);
        RDGeom::Point3D atom_pos = mol.getConformer().getAtomPos(atom->getIdx());

      }
    }
    rings.process_rings(mol);
    auto _rings_vec = rings.getRingsVector();

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
        Rings::find_center_and_normal(new_mol.getConformer(), ring, center, norm1, norm2);
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
    cleanup_predef(mol);

    if (result.has_unlisted_resnames) {
      clean_bonds(result.mol, result.mol.getConformer());
      perceive_bond_orders_obabel(result.mol);
      cleanup(result.mol);

      merge_bonds(mol, result.mol, result.atom_indices);
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
    st = read_structure_gz(file_name);
    gemmiStructureToRDKit(*mol, st, *conformer, false);
    mol->updatePropertyCache(false);
    mol->addConformer(conformer, true);
  }

public:
  explicit Luni(std::string file_name) : _cutoff(BONDED_NS_CUTOFF) {
    process_file(file_name);
    create_topology();
  }

  void create_topology() {
    const auto &conf = get_conformer();
    grid = FastNS(conf.getPositions(), _cutoff);
    bonded_nps = grid.self_search();

    Topology::compute_bonds(*mol, bonded_nps);
    topology.assign_atom_types(*mol);

    // topology.rings = std::move(top.rings);
    // topology = top;
  }

  // RDKit::RWMol &get_molecule() { return source->get_molecule(); }
  // const RDKit::RWMol &get_molecule() const { return source->get_molecule(); }

  const std::vector<RDGeom::Point3D> &positions(int confId = -1) const {
    return get_conformer(confId).getPositions();
  }

  const RDKit::Conformer &get_conformer(int id = -1) const {
    return mol->getConformer(id);
  }

  RDKit::RWMol &get_molecule() { return *mol; }
  const RDKit::RWMol &get_molecule() const { return *mol; }

  const double get_cutoff() const { return _cutoff; }

  const std::vector<AtomType> &get_atom_types() const { return topology.atom_types; }
  const RingDataVec &get_rings() const { return topology.rings_vec; }

  NSResults _find_neighbors(double cutoff = BONDED_NS_CUTOFF) {

    if (cutoff == _cutoff) {
      return bonded_nps;
    } else if (cutoff < _cutoff) {
      return bonded_nps.filter(cutoff);
    }

    grid.update_cutoff(cutoff);
    auto ns = grid.self_search();
    ns.m_luni = this;
    return ns;
  }

  template <typename T>
  Neighbors<T> find_neighbors(double cutoff) {

    if (fabs(cutoff - _cutoff) < 1e-6) {
      return {*this, bonded_nps.get_pairs(), bonded_nps.get_distances(), false};
    } else if (cutoff < _cutoff) {
      NSResults val = bonded_nps.filter(cutoff);
      return {*this, val.get_pairs(), val.get_distances(), false};
    }

    // const auto &conf = get_conformer();
    // auto _grid = FastNS(conf.getPositions(), cutoff);
    // NSResults new_ns = _grid.self_search();
    //
    // auto p = new_ns.get_pairs();
    // auto d = new_ns.get_distances();
    // std::vector<T> _p1;
    // for (size_t i = 0; i < p.size(); ++i) {
    //   _p1.push_back({p[i].first, p[i].second, d[i]});
    // }
    // Neighbors<T> _n1(*this, _p1);
    // return _n1;

    // FIX: update_cutoff is not working:
    grid.update_cutoff(cutoff);
    NSResults ns = grid.self_search();

    return Neighbors<T>(*this, ns.get_pairs(), ns.get_distances(), false);
  }

  std::vector<RDKit::MatchVectType>
  match_smarts_string(std::string sm, std::string atype = "",
                      bool log_values = false) const {

    // initialize ringinfo
    if (!mol->getRingInfo()->isInitialized()) {
      RDKit::MolOps::symmetrizeSSSR(*mol);
    }
    std::vector<RDKit::MatchVectType> match_list;
    auto sm_mol = RDKit::SmartsToMol(sm);
    mol->updatePropertyCache(false);
    RDKit::SubstructMatch(*mol, *sm_mol, match_list);
    return match_list;
  };

  // FIX: move to source
  const RDKit::Atom &get_atom(int index) const {
    return *mol->getAtomWithIdx(index);
  }

  const auto atoms() const { return mol->atoms(); }

  //! Returns the number of atoms in the molecule.
  const auto n_atoms() const { return mol->getNumAtoms(); }

  const auto coordinates() const {
    auto coords = mol->getConformer().getPositions();
    auto ccoords = reinterpret_cast<std::vector<std::vector<float>> &>(coords);
    return ccoords;
    // return reinterpret_cast<std::vector<std::vector<float>>&>(coords);
  }

  // manually create std::vector<std::vector<float>> from RDGeom::Point3D
  const auto coordinates2() const {
    auto coords = mol->getConformer().getPositions();
    std::vector<std::vector<double>> ccoords;
    for (const auto &coord : coords) {
      ccoords.push_back({coord.x, coord.y, coord.z});
    }
    return ccoords;
  }

  //! Returns the names of the atoms.
  const auto names() const {
    return atom_attrs_ref<std::string>(
        [](const RDKit::Atom *atom) -> const std::string & {
          auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(
              atom->getMonomerInfo());
          return info->getName();
        });
  }
  const auto _names() const {
    std::vector<std::string> names = std::vector<std::string>();
    for (const auto atom : mol->atoms()) {
      auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(
          atom->getMonomerInfo());
      names.push_back((std::string)info->getName());
    }
    return names;
  }

  const auto symbols() const {
    return atom_attrs<std::string>(
        [](const RDKit::Atom *atom) { return atom->getSymbol(); });
  }
  // const auto symbols() const {
  //   std::vector<std::string> symbols;
  //   for (const auto atom : source->get_molecule().atoms()) {
  //     symbols.push_back(atom->getSymbol());
  //   }
  //   return symbols;
  // }

  //! Returns the residue indices of the atoms.
  const auto indices() const {
    return atom_attrs<int>(
        [](const RDKit::Atom *atom) { return atom->getIdx(); });
  }

  //! Returns the atomic numbers of the atoms.
  const auto atomic_numbers() const {
    return atom_attrs<int>(
        [](const RDKit::Atom *atom) { return atom->getAtomicNum(); });
  }

  //! Returns the elements of the atoms.
  const auto elements() const {
    const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
    return atom_attrs<std::string>([&tbl](const RDKit::Atom *atom) {
      return tbl->getElementSymbol(atom->getAtomicNum());
    });
  }

  //! Returns the residue names of the atoms.
  const auto resnames() const {
    return atom_attrs_ref<std::string>(
        [](const RDKit::Atom *atom) -> const std::string & {
          auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(
              atom->getMonomerInfo());
          return info->getResidueName();
        });
  }

  const auto resids() const {
    return atom_attrs<int>([](const RDKit::Atom *atom) -> int {
      auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(
          atom->getMonomerInfo());
      return info->getResidueNumber();
    });
  }

  const auto resindices() const {
    return atom_attrs<int>([](const RDKit::Atom *atom) -> int {
      auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(
          atom->getMonomerInfo());
      return info->getSegmentNumber();
    });
  }

  const auto chainlabels() const {
    return atom_attrs_ref<std::string>(
        [](const RDKit::Atom *atom) -> const std::string & {
          auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(
              atom->getMonomerInfo());
          return info->getChainId();
        });
  }

  friend class NSResults;
  // friend Neighbors<AtomAtomPair>;
  // friend Neighbors<AtomRingPair>;

  template <typename T>
  friend class Neighbors;


private:
  template <typename T>
  std::vector<T> atom_attrs(std::function<T(const RDKit::Atom *)> func) const {
    std::vector<T> attrs;
    // auto &mol = source->get_molecule();
    attrs.reserve(mol->getNumAtoms());
    for (const auto atom : mol->atoms()) {
      attrs.push_back(func(atom));
    }
    return attrs;
  }

  template <typename T>
  std::vector<std::reference_wrapper<const T>>
  atom_attrs_ref(std::function<const T &(const RDKit::Atom *)> func) const {
    std::vector<std::reference_wrapper<const T>> attributes;
    // auto &mol = source->get_molecule();
    attributes.reserve(mol->getNumAtoms());
    for (const auto atom : mol->atoms()) {
      attributes.push_back(std::cref(func(atom)));
    }
    return attributes;
  }
};

} // namespace lahuta

#endif // LAHUTA_HPP
