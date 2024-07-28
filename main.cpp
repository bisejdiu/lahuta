#include "lahuta.hpp"
#include <gemmi/mmread_gz.hpp> // for read_structure_gz

#include <chrono>
#include <iostream>

#define TO_MS(d) std::chrono::duration_cast<std::chrono::milliseconds>(d)
#define T() std::chrono::high_resolution_clock::now()

using namespace gemmi;
using namespace RDKit;

int main(int argc, char const *argv[]) {
  auto startTotTime = std::chrono::high_resolution_clock::now();
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <cif file>" << std::endl;
    return 1;
  }
  std::string file_name = argv[1];

  auto load_start = T();
  Structure st = read_structure_gz(file_name);
  auto loadTime = TO_MS(T() - load_start).count();
  std::cout << "Load: " << loadTime << "ms" << std::endl;

  // Current API
  auto source = Lahuta::GemmiSource();
  source.process(st);

  Lahuta::Luni luni(source);

  // FIX: replace getNeighborPairsSize with just Size()
  std::cout << "Neighbors: " << luni.getNeighborResults().getNeighborPairsSize()
            << std::endl;

  auto r = luni.getNeighborResults();
  auto neighbors = r.getNeighbors();

  // float new_cutoff = 3.0;
  // r = luni.getNeighborResults().filterByDistance(new_cutoff);
  // std::cout << "Updated cutoff NP size: " << r.getNeighborPairsSize()
  //           << std::endl;
  //
  // float big_cutoff = 10.0;
  // r = luni.getNeighborResults().filterByDistance(big_cutoff);
  // std::cout << "filter (should be zero?): " << r.getNeighborPairsSize() << std::endl;
  //
  // r = luni.findNeighbors(big_cutoff);
  // std::cout << "Big cutoff NP size recomp: " << r.getNeighborPairsSize() << std::endl;
  // neighbors = r.getNeighbors();
  // for (size_t i = 0; i <= 10; ++i) {
  //   std::cout << neighbors[i].first << " " << neighbors[i].second << " "
  //             << r.distances[i] << std::endl;
  // }

  RDKit::RWMol *mol = &luni.getMolecule();
  int o1 = 0;
  int o2 = 0;
  for (auto bondIt = mol->beginBonds(); bondIt != mol->endBonds(); ++bondIt) {
    RDKit::Bond *bond = *bondIt;
    if (bond->getBondType() == RDKit::Bond::BondType::SINGLE) {
      o1++;
    } else if (bond->getBondType() == RDKit::Bond::BondType::DOUBLE) {
      o2++;
    }
  }
  std::cout << "1: " << o1 << " 2: " << o2 << " t: " << o1 + o2 << std::endl;
  std::cout << "Nr. Bonds: " << mol->getNumBonds() << std::endl;

  auto totTime = TO_MS(T() - startTotTime).count();
  std::cout << "Total Time: " << totTime << "ms" << std::endl;

  return 0;
}
