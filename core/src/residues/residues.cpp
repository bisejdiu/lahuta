#include <rdkit/GraphMol/MonomerInfo.h>

#include "find_rings.hpp"
#include "residues/definitions.hpp"
#include "residues/residues.hpp"

// clang-format off
namespace lahuta {

struct ResidueKeyHash {
  std::size_t operator()(const std::tuple<std::string, int, std::string, std::string> &t) const {
    std::size_t seed = 0;
    auto hash_combine = [&seed](const auto &val) {
      seed ^= std::hash<std::decay_t<decltype(val)>>{}(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    };

    hash_combine(std::get<0>(t));
    hash_combine(std::get<1>(t));
    hash_combine(std::get<2>(t));
    hash_combine(std::get<3>(t));

    return seed;
  }
};

bool Residues::build() {
  bool success = false;
  build_residues(mol_, success);
  return success;
}

Residues Residues::filter(std::function<bool(const Residue &)> predicate) const {
  Residues result(mol_);
  for (const auto &residue : residues_) {
    if (predicate(residue)) {
      result.residues_.push_back(residue);
    }
  }
  return result;
}

Residues Residues::filter(std::function<bool(const std::string &)> predicate) const {
  Residues result(mol_);
  for (const auto &residue : residues_) {
    if (predicate(residue.name)) {
      result.residues_.push_back(residue);
    }
  }
  return result;
}

template <typename ResultType>
std::vector<ResultType> Residues::map(std::function<ResultType(const Residue &)> func) const {
  std::vector<ResultType> result;
  result.reserve(residues_.size());
  for (const auto &residue : residues_) {
    result.push_back(func(residue));
  }
  return result;
}

void Residues::build_residues(const RDKit::RWMol &mol, bool &status) {
  std::unordered_map<std::tuple<std::string, int, std::string, std::string>, size_t, ResidueKeyHash> residue_index_map;

  std::vector<Residue> residues;
  for (const auto &atom : mol.atoms()) {
    // if (atom->getAtomicNum() == 1) continue;

    auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    if (!info) continue;

    int res_num          = info->getResidueNumber();
    std::string res_name = info->getResidueName();
    std::string chain_id = info->getChainId();

    // TODO: add a flag instead of the altloc
    /*std::string alt_loc = info->getAltLoc();*/
    /*bool has_alt_loc = !alt_loc.empty();*/

    auto key = std::make_tuple(chain_id, res_num, res_name, "");

    auto it = residue_index_map.find(key);
    if (it == residue_index_map.end()) {
      residues.emplace_back(chain_id, res_num, res_name, "");
      residue_index_map[key] = residues.size() - 1;
    }
    residues[residue_index_map[key]].atoms.push_back(atom);
  }

  for (std::size_t i = 0; i < residues.size(); ++i) {
    residues[i].idx = static_cast<unsigned>(i);
  }

  // residues are in the order how they are inserted into the map
  residues_ = std::move(residues);
  status = true;
}

std::vector<std::vector<int>> find_and_process_aromatic_residues(const RDKit::RWMol &mol, const Residues &residues) {

  using RingSize = definitions::arom_rings::RingSize;
  constexpr auto &AromRingSize = definitions::arom_rings::AromaticResiduesRingSizes;

  std::vector<std::vector<int>> ring_vector;
  for (const auto &residue : residues) {

    RingSize ring_size;
    if      (residue.name == "PHE") { ring_size = RingSize::RS_6; }
    else if (residue.name == "TYR") { ring_size = RingSize::RS_6; }
    else if (residue.name == "HIS") { ring_size = RingSize::RS_5; }
    else if (residue.name == "TRP") { ring_size = static_cast<RingSize>(RingSize::RS_5 | RingSize::RS_6); } 
    else {
      // linear search through remaining predefined residues
      auto it = std::find_if(AromRingSize.begin(), AromRingSize.end(), [&](const auto &e) { return e.first == residue.name; });
      if (it == std::end(AromRingSize)) continue;

      ring_size = it->second;
    }

    std::vector<int> ring_sizes = definitions::arom_rings::get_ringsizes(ring_size);
    for (int ring_size : ring_sizes) {
      if (residue.atoms.size() < static_cast<size_t>(ring_size)) continue;

      std::vector<int> processed_data = FastRingFinder::find_ring_in_residue(residues.get_mol(), residue, ring_size);
      if (processed_data.size() != ring_size) continue;

      ring_vector.push_back(std::move(processed_data));
    }
  }

  return ring_vector;
}

} // namespace lahuta
