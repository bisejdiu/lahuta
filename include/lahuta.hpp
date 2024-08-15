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

    if (!result.has_unlisted_resnames) {
      mol.updatePropertyCache(false);
      auto start = std::chrono::high_resolution_clock::now();
      cleanup_predef(mol);
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed = end - start;
      std::cout << "Cleanup predef: " << elapsed.count() << "s" << std::endl;

      // std::cout << "starting atom typing.. " << std::endl;
      auto start2 = std::chrono::high_resolution_clock::now();
      for (auto atom : mol.atoms()) {
        auto *info =
            static_cast<RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
        AtomType atom_type = get_atom_type(atom);
        std::cout << info->getResidueName() << " " << info->getName() << " "
        << atom_type_to_string(atom_type) << std::endl;
      }
      auto end2 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed2 = end2 - start2;
      std::cout << "Atom typing: " << elapsed2.count() << "s" << std::endl;

      return;
    }

    result.mol.updatePropertyCache(false);
    std::cout << "Unlisted residues: " << result.atom_indices.size()
              << std::endl;
    clean_bonds(result.mol, result.mol.getConformer());
    result.mol.updatePropertyCache(false);
    perceive_bond_orders_obabel(result.mol);
    result.mol.updatePropertyCache(false);
    cleanup(result.mol);

    mol.updatePropertyCache(false);
    cleanup_predef(mol);
    merge_bonds(mol, result.mol, result.atom_indices);

    mol.updatePropertyCache(false);

    // auto start2 = std::chrono::high_resolution_clock::now();
    // for (auto atom : mol.atoms()) {
    //   auto *info =
    //       static_cast<RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    //   AtomType atom_type = get_atom_type(atom);
    //   std::cout << info->getResidueName() << " " << info->getName() << " " <<
    //   atomTypeToString(atom_type) << std::endl;
    // }
    // auto end2 = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed2 = end2 - start2;
    // std::cout << "Atom typing: " << elapsed2.count() << "s" << std::endl;
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
