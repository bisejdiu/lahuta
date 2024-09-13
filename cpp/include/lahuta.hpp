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
#include "ob/clean_mol.hpp"
#include "ob/kekulize.h"

#define LAHUTA_VERSION "0.10.0"

namespace lahuta {

static float BONDED_NS_CUTOFF = 4.5;

// class ISource {
// public:
//   // FIX: RDKit provides ROMol and RWMol
//   virtual RDKit::RWMol &get_molecule() = 0;
//   virtual const RDKit::RWMol &get_molecule() const = 0;
//   virtual RDKit::Conformer &get_conformer(int id = -1) = 0;
//   virtual Structure &get_structure() = 0; // FIX: delete after test
//   virtual ~ISource() = default;
//
//   virtual void process(std::string file_name) = 0;
// };

class GemmiSource {
private:
  std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();
  Structure st;

public:
  explicit GemmiSource() = default;

  GemmiSource(const GemmiSource &source) {
    mol = std::make_unique<RDKit::RWMol>(*source.mol);
  }

  void process(std::string file_name) {
    RDKit::Conformer *conformer = new RDKit::Conformer();
    st = read_structure_gz(file_name);
    gemmiStructureToRDKit(*mol, st, *conformer, false);
    mol->updatePropertyCache(false);
    mol->addConformer(conformer, true);
  }

  RDKit::RWMol &get_molecule() { return *mol; }
  const RDKit::RWMol &get_molecule() const { return *mol; }

  RDKit::Conformer &get_conformer(int id = -1) {
    return mol->getConformer(id);
  }

  Structure &get_structure() { return st; }
};

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

struct RingData {
    std::vector<int> atom_ids;
    RDGeom::Point3D center;
    RDGeom::Point3D norm1;
    RDGeom::Point3D norm2;

    RingData() = default;
    RingData(RDGeom::Point3D center, RDGeom::Point3D norm1, RDGeom::Point3D norm2,
             std::vector<int> atom_ids)
        : center(center), norm1(norm1), norm2(norm2), atom_ids(atom_ids) {}
};

struct RingDataVec {
  std::vector<RingData> rings;

  RingDataVec() = default;
  RingDataVec(RingDataVec &&other) noexcept : rings(std::move(other.rings)) {}
  RingDataVec &operator=(RingDataVec &&other) noexcept {
    if (this != &other) {
      rings = std::move(other.rings);
    }
    return *this;
  }

  // RingDataVec(std::vector<RingData> rings) : rings(std::move(rings)) {}
  // // copy constructor
  RingDataVec(const RingDataVec &other) : rings(other.rings) {}
  // copy = operator
  RingDataVec &operator=(const RingDataVec &other) {
    if (this != &other) {
      rings = other.rings;
    }
    return *this;
  }
  
  std::vector<std::vector<double>> centers() const {
    std::vector<std::vector<double>> centers;
    for (const auto &ring : rings) {
      centers.push_back({ring.center.x, ring.center.y, ring.center.z});
    }
    return centers;
  }

  RDGeom::POINT3D_VECT centers_rkdit() const {
    RDGeom::POINT3D_VECT centers;
    for (const auto &ring : rings) {
      centers.push_back(ring.center);
    }
    return centers;
  }

  std::vector<std::vector<double>> norm1() const {
    std::vector<std::vector<double>> normals;
    for (const auto &ring : rings) {
      normals.push_back({ring.norm1.x, ring.norm1.y, ring.norm1.z});
    }
    return normals;
  }

  std::vector<std::vector<double>> norm2() const {
    std::vector<std::vector<double>> normals;
    for (const auto &ring : rings) {
      normals.push_back({ring.norm2.x, ring.norm2.y, ring.norm2.z});
    }
    return normals;
  }
};

class Rings {
public:
  struct ResId {
    std::string chain;
    int res_num;
    std::string res_name;

    bool operator==(const ResId &other) const {
      return chain == other.chain && res_num == other.res_num &&
             res_name == other.res_name;
    }
  };

  struct ResIdHash {
    std::size_t operator()(const ResId &id) const {
      return std::hash<std::string>()(id.chain) ^
             (std::hash<int>()(id.res_num) << 1) ^
             (std::hash<std::string>()(id.res_name) << 2);
    }
  };


  using RingMap = std::unordered_map<ResId, RingData, ResIdHash>;
  using AtomOrderMap =
      std::unordered_map<std::string, std::vector<std::string>>;

  Rings() { initializeAtomOrders(); }

  void add_ring_atom(const RDKit::Atom *atom) {
    auto *info =
        static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    if (!info)
      return;

    ResId res_id{info->getChainId(), info->getResidueNumber(),
                 info->getResidueName()};
    const bool is_trp = info->getResidueName() == "TRP";
    if (is_trp) {
      res_id.res_name = "TRP5";
      rings_[res_id].atom_ids.push_back(atom->getIdx());
      res_id.res_name = "TRP6";
    }
    rings_[res_id].atom_ids.push_back(atom->getIdx());
  }

  void process_rings(const RDKit::ROMol &mol) {
    const RDKit::Conformer &conf = mol.getConformer();
    for (auto &[res_id, ring_data] : rings_) {
      if (!ring_data.atom_ids.empty()) {
        order_atoms(res_id.res_name, ring_data.atom_ids, mol);
        find_center_and_normal(conf, ring_data.atom_ids, ring_data.center,
                               ring_data.norm1, ring_data.norm2);
      }
    }
  }

  RingData getRing(const ResId &res_id) { return rings_[res_id]; }

  const RingDataVec getRingsVector() {
    RingDataVec ring_data;
    for (auto &[res_id, data] : rings_) {
      ring_data.rings.push_back(data);
    }
    return ring_data;
  }


  const RingMap &getRings() const { return rings_; }

private:
  RingMap rings_;
  AtomOrderMap atom_orders_;

  void initializeAtomOrders() {
    atom_orders_["PRO"] = {"CG", "CB", "CA", "N", "CD"};
    atom_orders_["TYR"] = {"CE2", "CZ", "CE1", "CD1", "CG", "CD2"};
    atom_orders_["HIS"] = {"CD2", "NE2", "CE1", "ND1", "CG"};
    atom_orders_["PHE"] = {"CE1", "CZ", "CE2", "CD2", "CG", "CD1"};
    atom_orders_["TRP5"] = {"NE1", "CE2", "CD2", "CG", "CD1"};
    atom_orders_["TRP6"] = {"CZ2", "CH2", "CZ3", "CE3", "CD2", "CE2"};
  }

  void order_atoms(const std::string &res_name, std::vector<int> &atom_ids,
                   const RDKit::ROMol &mol) {
    if (atom_orders_.find(res_name) == atom_orders_.end())
      return;

    const auto &expected_order = atom_orders_[res_name];
    std::vector<int> ordered_ids;
    ordered_ids.reserve(expected_order.size());

    for (const auto &exptected_name : expected_order) {
      auto it = std::find_if(atom_ids.begin(), atom_ids.end(), [&](int idx) {
        const RDKit::Atom *atom = mol.getAtomWithIdx(idx);
        auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(
            atom->getMonomerInfo());
        return info && info->getName() == exptected_name;
      });
      if (it != atom_ids.end()) {
        ordered_ids.push_back(*it);
      }
    }
    atom_ids = std::move(ordered_ids);
    // NOTE: Data can have missing atoms, which is currently not considered.
    // This probably only matters for PRO since it's not aromatic and hence not
    // planar.
  }

public:
  static void find_center_and_normal(const RDKit::Conformer &conf,
                                     const std::vector<int> &ring,
                                     RDGeom::Point3D &center,
                                     RDGeom::Point3D &norm1,
                                     RDGeom::Point3D &norm2) {

    center = std::accumulate(ring.begin(), ring.end(), center,
                             [&conf](const RDGeom::Point3D &sum, int idx) {
                               return sum + conf.getAtomPos(idx);
                             });
    center /= static_cast<double>(ring.size());

    for (size_t i = 0; i < ring.size(); ++i) {
      RDGeom::Point3D v1 = conf.getAtomPos(ring[i]) - center;
      RDGeom::Point3D v2 =
          conf.getAtomPos(ring[(i + 1) % ring.size()]) - center;
      norm1 += v1.crossProduct(v2);
    }
    norm1 /= static_cast<double>(ring.size());
    norm1.normalize();
    norm2 = -norm1;
  }
};

// struct Topology {
//   std::vector<AtomType> atom_types;
//   std::vector<std::unique_ptr<RingData>> rings;
// };


class Topology {
public:
  std::vector<AtomType> atom_types;
  RingDataVec rings_vec;

public:
  static void cleanup_predef(RDKit::RWMol &mol) { basic_mol_cleanup(mol); }

  static void cleanup(RDKit::RWMol &mol) {
    basic_mol_cleanup(mol);

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
    // assign_atom_types(mol);
  }
};




class Luni {
private:
  // std::shared_ptr<GemmiSource> source;
  std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();
  Structure st;
  NSResults bonded_nps;
  float _cutoff;
  FastNS grid;
  Topology topology;
  // std::vector<AtomType> atom_types;

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

    // atom_types = topology.get_atom_types();


    // topology.rings = std::move(top.rings);
    // topology = top;

    // debug atom 608
    // auto atom = source->get_molecule().getAtomWithIdx(608);
    // auto *info = static_cast<RDKit::AtomPDBResidueInfo
    // *>(atom->getMonomerInfo()); auto atom_type = atom_types[608]; std::cout
    // << "Atom 608: " << info->getResidueName() << " " << info->getName() << "
    // " << atom_type_to_string(atom_type) << std::endl; debug all ASP residues:
    // for (auto &atom: source->get_molecule().atoms()) {
    //   auto *info = static_cast<RDKit::AtomPDBResidueInfo
    //   *>(atom->getMonomerInfo()); if (info->getResidueName() == "ASP") {
    //     auto atom_type = atom_types[atom->getIdx()];
    //     std::cout << "ASP: " << info->getResidueName() << " " <<
    //     info->getName() << " " << atom_type_to_string(atom_type) <<
    //     std::endl;
    //   }
    // }
  }

  // RDKit::RWMol &get_molecule() { return source->get_molecule(); }
  // const RDKit::RWMol &get_molecule() const { return source->get_molecule(); }

  const std::vector<RDGeom::Point3D> &positions(int confId = -1) const {
    return get_conformer(confId).getPositions();
  }

  const RDKit::Conformer &get_conformer(int id = -1) const {
    return mol->getConformer(id);
  }
  // RDKit::Conformer &get_conformer(int id = -1) {
  //   return mol->getConformer(id);
  // }
  RDKit::RWMol &get_molecule() { return *mol; }
  const RDKit::RWMol &get_molecule() const { return *mol; }

  const double get_cutoff() const { return _cutoff; }

  const std::vector<AtomType> &get_atom_types() const { return topology.atom_types; }
  const RingDataVec &get_rings() const { return topology.rings_vec; }

  NSResults find_neighbors(double cutoff = BONDED_NS_CUTOFF) {

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
