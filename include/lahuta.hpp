#include <gemmi/mmread_gz.hpp> // for read_structure_gz
#include <rdkit/GraphMol/BondIterators.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>

#include "bond_order.hpp"
#include "bonds.hpp"
#include "convert.hpp"
#include "gemmi/model.hpp"
#include "nsgrid.hpp"
#include "ob/clean_mol.hpp"
#include "ob/kekulize.h"

#include "atom_types.hpp"

#include <chrono>
#include <memory>
#include <optional>
#include <string>

#define LAHUTA_VERSION "0.10.0"

namespace Lahuta {


inline std::vector<AtomType> AtomTypeMatch(RDKit::ROMol &mol, SubstructMatchParameters &params) {
  static std::array<RDKit::ROMol *, std::size(AtomTypeSMARTS)> patterns = [] {
    std::array<RDKit::ROMol *, std::size(AtomTypeSMARTS)> temp{};
    for (size_t i = 0; i < std::size(AtomTypeSMARTS); ++i) {
      temp[i] = RDKit::SmartsToMol(AtomTypeSMARTS[i].first);
    }
    return temp;
  }();

  std::vector<AtomType> types = {mol.getNumAtoms(), AtomType::NONE};
  // types.reserve(mol.getNumAtoms());
  for (size_t i = 0; i < std::size(AtomTypeSMARTS); ++i) {
    const auto &[smarts, atom_type] = AtomTypeSMARTS[i];
    RDKit::ROMol *pattern = patterns[i];

    SubStrMatches matchList;
    RDKit::SubstructMatch(mol, *pattern, matchList);

    for (const auto &match : matchList) {
      // std::cout << "Match: " << match[0].second << " " << atom_type_to_string(atom_type) << std::endl;
      types[match[0].second] |= atom_type;
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

    // auto *i = static_cast<RDKit::AtomPDBResidueInfo*>(at->getMonomerInfo());
    mol.updatePropertyCache(false);
    cleanup_predef(mol);



    std::vector<int> invalid_indices; 
    invalid_indices.reserve(mol.getNumAtoms());
    for (auto atom : mol.atoms()) {
      AtomType atom_type = get_atom_type(atom);
      if (atom_type == AtomType::INVALID && atom->getAtomicNum() != 1) {
        invalid_indices.push_back(atom->getIdx());
      }
    }

    std::cout << "sizes: " << result.atom_indices.size() << " " << invalid_indices.size() << std::endl;
    auto with_bonds = true;
    auto new_mol = rdMolFromRDKitMol(mol, invalid_indices, with_bonds);
    std::cout << "created new mol\n";
    SubstructMatchParameters params;
    params.maxMatches = new_mol.getNumAtoms();
    std::cout << "sizes all: " << result.atom_indices.size() << " " << invalid_indices.size() << " " << new_mol.getNumAtoms() << std::endl;
    
    // initialize ringinfo
    if (!new_mol.getRingInfo()->isInitialized()) {
      RDKit::MolOps::symmetrizeSSSR(new_mol);
    }
    auto vec = AtomTypeMatch(new_mol, params);
    for (auto i = 0; i < vec.size(); ++i) {
      // AtomType atom_type = get_atom_type(atom);
      AtomType atom_type = vec[i];
      auto atom = new_mol.getAtomWithIdx(i);
      if (atom_type != AtomType::INVALID) {
        std::cout << "x atom: " << atom->getIdx() << " " << atom->getSymbol() << " " << atom_type_to_string(vec[i]) << std::endl;
        // std::cout << "x atom: " << atom_type_to_string(atom_type) << std::endl;
      }
    }

    if (result.has_unlisted_resnames) {
      std::cout << "Unlisted residues: " << result.atom_indices.size()
                << std::endl;
      clean_bonds(result.mol, result.mol.getConformer());
      perceive_bond_orders_obabel(result.mol);
      cleanup(result.mol);
 
      // FIX: test if insertMol is a solution. In our case we are not so much merging or inserting, but 
      // rather adding bonds on mol based on result.mol connectivity and bond types. 
      merge_bonds(mol, result.mol, result.atom_indices);
    }
  }
};

class Luni {
private:
  std::unique_ptr<ISource> source;
  NSResults neighbors;
  float _cutoff;
  FastNS grid;

public:
  explicit Luni(std::string file_name) : _cutoff(4.5) {
    source = std::make_unique<GemmiSource>();
    source->process(file_name);
    init_and_compute_bonds();
  }

  void init_and_compute_bonds() {
    const auto &conf = source->get_conformer();
    grid = FastNS(conf.getPositions(), _cutoff);
    neighbors = grid.self_search();
    BondComputation::compute_bonds(source->get_molecule(),
                                   source->get_structure(), neighbors);
  }

  RDKit::RWMol &get_molecule() { return source->get_molecule(); }
  const RDKit::RWMol &get_molecule() const { return source->get_molecule(); }
  const std::vector<RDGeom::Point3D> &positions(int confId = -1) const {
    return source->get_conformer(confId).getPositions();
  }
  const NeighborPairs &get_neighbors() const {
    return neighbors.get_neighbors();
  }
  const std::vector<float> &get_distances() const {
    return neighbors.get_distances();
  }
  const double get_cutoff() const { return _cutoff; }

  NSResults find_neighbors(std::optional<double> cutoff) {
    auto value = cutoff.value_or(_cutoff);

    if (value == _cutoff) {
      return neighbors;
    } else if (value < _cutoff) {
      return neighbors.filter(value);
    }

    grid.update_cutoff(value);
    return grid.self_search();
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
  }
};

} // namespace Lahuta
