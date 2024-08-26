#ifndef LAHUTA_HPP
#define LAHUTA_HPP

#include <chrono>
#include <memory>
#include <optional>
#include <string>
#include <array>
#include <vector>

#include <gemmi/mmread_gz.hpp> // for read_structure_gz
#include <rdkit/GraphMol/BondIterators.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>

#include <boost/range/iterator_range.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "bond_order.hpp"
#include "bonds.hpp"
#include "convert.hpp"
#include "gemmi/model.hpp"
#include "nsgrid.hpp"
#include "ob/clean_mol.hpp"
#include "ob/kekulize.h"

// #include "atom_types.hpp"


#define LAHUTA_VERSION "0.10.0"

namespace lahuta {

inline std::vector<AtomType> match_atom_types(RDKit::ROMol &mol) {
  static std::array<RDKit::ROMol *, std::size(AtomTypeSMARTS)> patterns = [] {
    std::array<RDKit::ROMol *, std::size(AtomTypeSMARTS)> temp{};
    for (size_t i = 0; i < std::size(AtomTypeSMARTS); ++i) {
      temp[i] = RDKit::SmartsToMol(AtomTypeSMARTS[i].first);
    }
    return temp;
  }();

  SubstructMatchParameters params;
  params.maxMatches = mol.getNumAtoms();

  std::vector<AtomType> types = {mol.getNumAtoms(), AtomType::NONE};
  for (size_t i = 0; i < std::size(AtomTypeSMARTS); ++i) {
    const auto &[smarts, atom_type] = AtomTypeSMARTS[i];
    RDKit::ROMol *pattern = patterns[i];

    SubStrMatches matchList;
    RDKit::SubstructMatch(mol, *pattern, matchList);

    for (const auto &match : matchList) {
      for (const auto &pair : match) {
        types[pair.second] |= atom_type;
      }
      // auto *atom = mol.getAtomWithIdx(match[0].second);
        // if (atom->getAtomicNum() == 26) {
        //   auto *info = static_cast<RDKit::AtomPDBResidueInfo*>(atom->getMonomerInfo());
        //   std::cout << "-> Fe: " << info->getResidueName() << " "
        //             << info->getName() << " " << atom_type_to_string(atom_type)
        //             << std::endl;
        // }
      // types[match[0].second] |= atom_type;
    }
  }

  return types;
}



class ISource {
public:
  // FIX: RDKit provides ROMol and RWMol
  virtual RDKit::RWMol &get_molecule() = 0;
  virtual const RDKit::RWMol &get_molecule() const = 0;
  virtual RDKit::Conformer &get_conformer(int id = -1) = 0;
  virtual Structure &get_structure() = 0; // FIX: delete after test
  virtual ~ISource() = default;

  virtual void process(std::string file_name) = 0;
};

class GemmiSource : public ISource {
private:
  std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();
  Structure st;

public:
  explicit GemmiSource() = default;

  GemmiSource(const GemmiSource &source) {
    mol = std::make_unique<RDKit::RWMol>(*source.mol);
  }

  void process(std::string file_name) override {
    RDKit::Conformer *conformer = new RDKit::Conformer();
    st = read_structure_gz(file_name);
    gemmiStructureToRDKit(*mol, st, *conformer, false);
    mol->updatePropertyCache(false);
    mol->addConformer(conformer, true);
  }

  RDKit::RWMol &get_molecule() override { return *mol; }
  const RDKit::RWMol &get_molecule() const override { return *mol; }

  RDKit::Conformer &get_conformer(int id = -1) override {
    return mol->getConformer(id);
  }

  Structure &get_structure() override { return st; }
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

class BondComputation {
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

  static std::vector<AtomType> assign_atom_types(RDKit::RWMol &mol) {

    std::vector<AtomType> indices(mol.getNumAtoms(), AtomType::NONE);
    std::vector<int> invalid_indices;

    // First pass: populate `indices` and track invalid non-hydrogen atoms.
    for (auto atom : mol.atoms()) {
      AtomType atom_type = get_atom_type(atom);
      indices[atom->getIdx()] = atom_type;
      if (atom_type == AtomType::INVALID && atom->getAtomicNum() != 1) {
        invalid_indices.push_back(atom->getIdx());
      } 
    }

    // Handle invalid indices using a new RDKit object.
    if (!invalid_indices.empty()) {
      std::cout << "Invalid atom types: " << invalid_indices.size() << std::endl;
      auto new_mol = filter_with_bonds(mol, invalid_indices);
      std::cout << "New mol atoms: " << new_mol.getNumAtoms() << std::endl;
      if (!mol.getRingInfo()->isInitialized()) {
        RDKit::MolOps::symmetrizeSSSR(new_mol);
      }
      auto vec = match_atom_types(new_mol);

      for (size_t i = 0; i < invalid_indices.size(); ++i) {
        indices[invalid_indices[i]] = vec[i];
      }
    }

    // for (size_t i = 0; i < indices.size(); ++i) {
    // for (size_t i = 0; i < 10; ++i) {
    //   auto atom = mol.getAtomWithIdx(i);
    //   auto *info =
    //       static_cast<RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    //
    //   if (atom->getAtomicNum() != 1) {
    //     std::cout << info->getResidueName() << " " << info->getName() << " "
    //               << atom_type_to_string(indices[i]) << std::endl;
    //   }
    // }

    return indices;
  }


  static void compute_bonds(RDKit::RWMol &mol, Structure &st,
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

    // assign_atom_types(mol);
  }
};



class Luni {
private:
  std::shared_ptr<ISource> source;
  NSResults bonded_neighbors;
  float _cutoff;
  FastNS grid;
  std::vector<AtomType> atom_types;

public:
  explicit Luni(std::string file_name) : _cutoff(4.5) {
    source = std::make_unique<GemmiSource>();
    source->process(file_name);
    init_and_compute_bonds();
  }

  void init_and_compute_bonds() {
    const auto &conf = source->get_conformer();
    grid = FastNS(conf.getPositions(), _cutoff);
    bonded_neighbors = grid.self_search();
    BondComputation::compute_bonds(source->get_molecule(),
                                   source->get_structure(), bonded_neighbors);
    atom_types = BondComputation::assign_atom_types(source->get_molecule());
  }

  RDKit::RWMol &get_molecule() { return source->get_molecule(); }
  const RDKit::RWMol &get_molecule() const { return source->get_molecule(); }
  const std::vector<RDGeom::Point3D> &positions(int confId = -1) const {
    return source->get_conformer(confId).getPositions();
  }
  auto &get_neighbors() const {
    return bonded_neighbors.get_neighbors();
  }
  // const std::vector<float> &get_distances() const {
  //   return neighbors.get_distances();
  // }
  const double get_cutoff() const { return _cutoff; }

  const std::vector<AtomType> &get_atom_types() const { return atom_types; }
  
  const RDKit::Atom &get_atom(int index) const {
    return *source->get_molecule().getAtomWithIdx(index);
  }

  NSResults filter_by_atom_type(AtomType type, int partner) {
    _NeighborPairs filtered;
    std::vector<float> distances;
    for (size_t i = 0; i < bonded_neighbors.size(); ++i) {
      if (partner == 0) {
        if (has(atom_types[bonded_neighbors.get_neighbors()[i].first], type)) {
          filtered.push_back(bonded_neighbors.get_neighbors()[i]);
          distances.push_back(bonded_neighbors.get_distances()[i]);
        }
      } else if (partner == 1) {
        if (has(atom_types[bonded_neighbors.get_neighbors()[i].second], type)) {
          filtered.push_back(bonded_neighbors.get_neighbors()[i]);
          distances.push_back(bonded_neighbors.get_distances()[i]);
        }
      } else {
        std::cerr << "Invalid partner: " << partner << std::endl;
      }
    }
    return NSResults(*this, std::move(filtered), std::move(distances));
  }

  NSResults find_neighbors(double cutoff = 4.5) {

    if (cutoff == _cutoff) {
      return bonded_neighbors;
    } else if (cutoff < _cutoff) {
      return bonded_neighbors.filter(cutoff);
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
    if (!source->get_molecule().getRingInfo()->isInitialized()) {
      RDKit::MolOps::symmetrizeSSSR(source->get_molecule());
    }
    std::vector<RDKit::MatchVectType> match_list;
    auto sm_mol = RDKit::SmartsToMol(sm);
    source->get_molecule().updatePropertyCache(false);
    RDKit::SubstructMatch(source->get_molecule(), *sm_mol, match_list);
    return match_list;
  };


  friend class NSResults;

};

struct Pairs {
    std::pair<int, int> neighbor_pair;
    float distance;
};

} // namespace Lahuta

#endif // LAHUTA_HPP
