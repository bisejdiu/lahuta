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
