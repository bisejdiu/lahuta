#include <rdkit/GraphMol/BondIterators.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <gemmi/mmread_gz.hpp> // for read_structure_gz

#include "bonds.hpp"
#include "bond_order.hpp"
#include "convert.hpp"
#include "nsgrid.hpp"
#include "ob/clean_mol.hpp"

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
  virtual ~ISource() = default;

  virtual void process(std::string file_name) = 0;
};

class GemmiSource : public ISource {
private:
  std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();

public:
  explicit GemmiSource() = default;

  GemmiSource(const GemmiSource &source) {
    mol = std::make_unique<RDKit::RWMol>(*source.mol);
  }

  void process(std::string file_name) override {
    RDKit::Conformer *conformer = new RDKit::Conformer();
    Structure st = read_structure_gz(file_name);
    gemmiStructureToRDKit(*mol, st, *conformer, false);
    mol->addConformer(conformer, true);
  }

  RDKit::RWMol &get_molecule() override { return *mol; }
  const RDKit::RWMol &get_molecule() const override { return *mol; }

  RDKit::Conformer &get_conformer(int id = -1) override {
    return mol->getConformer(id);
  }
};

class BondComputation {
public:
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

  static void compute_bonds(RDKit::RWMol &mol, const NSResults &neighborResults) {
    auto result = assign_bonds(mol, neighborResults);
    if (!result.has_unlisted_resnames) {
      return;
    }
    clean_bonds(result.mol, result.mol.getConformer());
    result.mol.updatePropertyCache(false);
    perceive_bond_orders_obabel(result.mol);
    merge_bonds(mol, result.mol, result.atom_indices);
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
    BondComputation::compute_bonds(source->get_molecule(), neighbors);
  }

  RDKit::RWMol &get_molecule() { return source->get_molecule(); }
  const RDKit::RWMol &get_molecule() const { return source->get_molecule(); }
  const std::vector<RDGeom::Point3D> &positions(int confId = -1) const {
    return source->get_conformer(confId).getPositions();
  }
  const NeighborPairs &get_neighbors() const { return neighbors.get_neighbors(); }
  const std::vector<float> &get_distances() const { return neighbors.get_distances(); }
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
  int match_smarts_string(std::string sm) const {
    std::vector<RDKit::MatchVectType> match_list;
    auto sm_mol = RDKit::SmartsToMol(sm);
    source->get_molecule().updatePropertyCache(false);
    // NOTE: do H make a difference
    // std::cout << "smart string processing: " << sm_mol->getNumAtoms() << std::endl;
    RDKit::SubstructMatch(source->get_molecule(), *sm_mol, match_list);
    return match_list.size();
    // std::vector<int> results(match_list.size());
    // for (auto &match : match_list) {
    //   auto atom = source->get_molecule().getAtomWithIdx(match[0].second);
    //   // std::cout << atom->getIdx() << 
    //   results.push_back(match[0].second);
    // }
    // return 0;
  }
};

} // namespace Lahuta
