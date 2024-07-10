#include <gemmi/mmread_gz.hpp> // for read_structure_gz

#include "bonds.hpp"  // for perceiveBonds
#include "conv.hpp"   // for gemmiStructureToRDKit
#include "nsgrid.hpp" // for FastNS

#include <chrono>

#include "bond_order.hpp"
#include "bond_utils.hpp" // for CleanUpMolecule

#include "bond_table/bonds.hpp"
#include "bond_table/table.hpp"

#include <iostream>

using namespace gemmi;
using namespace RDKit;

std::vector<BondInfo> _findBondsDeconstructed(Structure &st, Model &model,
                                              double maxRadius) {
  NeighborSearch ns(model, st.cell, maxRadius);
  ns.populate();

  std::vector<int> flags, order, key;
  std::vector<BondInfo> ret;

  initialize_bond_order_table();

  for (const Chain &chain : model.chains) {
    for (const Residue &res : chain.residues) {
      for (const gemmi::Atom &atom : res.atoms) {
        auto indices = model.get_indices(&chain, &res, &atom);
        Position pos = atom.pos;
        char altloc = atom.altloc;
        double thresholdA = getElementThreshold(atom.element.atomic_number());

        // NOTE: still do not understand the condition below
        // Could it be a consequence of the way the NeighborSearch is
        // implemented (e.g. no symmetry check for the pairs we get? )
        ns.for_each(
            pos, altloc, maxRadius,
            [&](const NeighborSearch::Mark &m, double dist_sq) {
              // Ensure that (1) we don't double count by
              // only considering atoms that are in the same residue or
              // in a residue that came before the current residue in
              // the chain and (2) we don't consider atoms that are too
              // close to each other (dist_sq < 1)
              if (indices[0] > m.chain_idx ||
                  (indices[0] == m.chain_idx &&
                   (indices[1] > m.residue_idx ||
                    (indices[1] == m.residue_idx && dist_sq < 1)))) {

                // print all conditions to debug
                // const_CRA cra1 = {&chain, &res, &atom};
                // const_CRA cra2 = m.to_cra(model);
                // printf("cra1 and cra2: %s %s\n", atom_str(cra1).c_str(),
                //        atom_str(cra2).c_str());
                // printf("indices[0] < m.chain_idx %d %d\n", indices[0],
                //        m.chain_idx);
                // printf("indices[1] < m.residue_idx %d %d\n", indices[1],
                //         m.residue_idx);
                // printf("dist_sq < 1 %f\n", dist_sq);
                //
                //
                // printf("Skipping %s %s %s %d %d %f\n", chain.name.c_str(),
                //        res.name.c_str(), atom.name.c_str(), m.chain_idx,
                //        m.residue_idx, dist_sq);
                return;
              }

              const_CRA cra1 = {&chain, &res, &atom};
              const_CRA cra2 = m.to_cra(model);

              // skip B-A bonds and only leave A-B bonds
              if (cra1.atom->serial >= cra2.atom->serial) {
                return;
              }

              // ret.push_back({cra1, cra2, m.image_idx, dist_sq});

              const gemmi::Atom *neighborAtom = cra2.atom;

              int neighborElemIdx = getElementIdx(neighborAtom->element.name());
              // neighborAtom.element.atomic_number();
              double thresholdB = getElementThreshold(neighborElemIdx);
              //
              double dist = std::sqrt(dist_sq);
              double pairingThreshold =
                  getPairingThreshold(atom.element.atomic_number(),
                                      neighborElemIdx, thresholdA, thresholdB);
              //
              if (dist <= pairingThreshold) {
                int order = get_intra_bond_order(res.name, atom.name,
                                                 neighborAtom->name);

                ret.push_back({cra1, cra2, m.image_idx, dist_sq, order});
              }
            });
      }
    }
  }
  return ret;
};
int main(int argc, char const *argv[]) {
  auto startTotTime = std::chrono::high_resolution_clock::now();
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <cif file>" << std::endl;
    return 1;
  }
  std::string file_name = argv[1];
  bool flag_b = false;

  for (int i = 2; i < argc; ++i) {
    if (std::string(argv[i]) == "-b") {
      flag_b = true;
    }
  }

  auto start = std::chrono::high_resolution_clock::now();
  Structure st = read_structure_gz(file_name);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Time to load structure: " << elapsed.count() * 1000 << " ms"
            << std::endl;

  // Model model = st.first_model();
  // // find_bonds(st, model);
  // auto ret = _findBondsDeconstructed(st, model, 4.0);
  // std::cout << "0. Number of bonds found: " << ret.size() << std::endl;
  // for (auto &bond : ret) {
  //   bond.print(st.cell);
  // }

  start = std::chrono::high_resolution_clock::now();
  RDKit::Conformer *conf = new RDKit::Conformer();
  RDKit::RWMol mol = gemmiStructureToRDKit(st, *conf, false);
  mol.addConformer(conf, true);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to convert to RDKit: " << elapsed.count() * 1000 << " ms"
            << std::endl;

  double cutoff = 4.5; // 01;
  std::vector<RDGeom::Point3D> atom_coords = conf->getPositions();
  std::cout << "Number of coords vs Number of atoms: " << atom_coords.size() << " " << mol.getNumAtoms() << std::endl;
  FastNS grid(atom_coords, cutoff);
  auto results = grid.selfSearch();
  std::cout << "Number of neighbors: " << results.getNeighborPairsSize()
            << std::endl;

  if (!flag_b) {
    // perceives bonds
    start = std::chrono::high_resolution_clock::now();
    perceiveBonds(mol, results, 0.45);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Time to perceive bonds: " << elapsed.count() * 1000 << " ms"
              << std::endl;
    std::cout << "Number of bonds found using NSGrid: " << mol.getNumBonds()
              << std::endl;

    // NOTE: should we clean up bonds before or after adding the struct conns?
    mol.updatePropertyCache(false);
    // Clean up incorrect bonds
    start = std::chrono::high_resolution_clock::now();
    CleanUpMolecule(mol, *conf);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Time to clean up molecule: " << elapsed.count() * 1000
              << "ms " << std::endl;

    for (gemmi::Connection &conn : st.connections) {
      // FIX: Need to iterate over all models
      gemmi::Atom *a1 = st.first_model().find_cra(conn.partner1).atom;
      gemmi::Atom *a2 = st.first_model().find_cra(conn.partner2).atom;

      if (mol.getBondBetweenAtoms(a1->serial - 1, a2->serial - 1) == nullptr) {
        mol.addBond((unsigned int)a1->serial - 1, (unsigned int)a2->serial - 1,
                    RDKit::Bond::BondType::SINGLE);
        continue;
      }
    }

    for (auto atomIt = mol.beginAtoms(); atomIt != mol.endAtoms(); ++atomIt) {
      auto atom = *atomIt;
      atom->setHybridization(HybridizationType::SP);
    }
    start = std::chrono::high_resolution_clock::now();
    PerceiveBondOrders(mol);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Time to perceive bond orders: " << elapsed.count() * 1000
              << " ms" << std::endl;
  }
  // for (gemmi::Connection &conn : st.connections) {
  //   // FIX: Need to iterate over all models
  //   gemmi::Atom *a1 = st.first_model().find_cra(conn.partner1).atom;
  //   gemmi::Atom *a2 = st.first_model().find_cra(conn.partner2).atom;
  //
  //   if (mol.getBondBetweenAtoms(a1->serial - 1, a2->serial - 1) == nullptr) {
  //     mol.addBond((unsigned int)a1->serial - 1, (unsigned int)a2->serial - 1,
  //                 RDKit::Bond::BondType::SINGLE);
  //     continue;
  //   }
  // }
  //
  // std::cout << "Number of bonds w/ _struct_conns " << mol.getNumBonds()
  //           << std::endl;

  // get ROMol from RWMol using cast
  // auto _mol = static_cast<RDKit::ROMol*>(new RDKit::RWMol(mol));

  if (flag_b) {
    start = std::chrono::high_resolution_clock::now();
    findBondsDeconstructedRDKit(mol, results);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Table lookup bonds: " << elapsed.count() * 1000 << " ms"
              << std::endl;
  }

  std::string bondOrders = "";
  for (auto bondIt = mol.beginBonds(); bondIt != mol.endBonds(); ++bondIt) {
    RDKit::Bond *bond = *bondIt;
    auto res1 = mol.getAtomWithIdx(bond->getBeginAtomIdx())->getMonomerInfo();
    auto res2 = mol.getAtomWithIdx(bond->getEndAtomIdx())->getMonomerInfo();

    AtomPDBResidueInfo *residue1 = dynamic_cast<AtomPDBResidueInfo *>(res1);
    AtomPDBResidueInfo *residue2 = dynamic_cast<AtomPDBResidueInfo *>(res2);

    bondOrders += std::to_string(bond->getBeginAtomIdx()) + " " +
                  // res1->getName() + " " + residue1->getResidueName() + " " + residue1->getAltLoc() + " " +
                  std::to_string(bond->getEndAtomIdx()) + " " + 
                  // res2->getName() + " " + residue2->getResidueName() + " " + residue2->getAltLoc() + " " +
                  std::to_string(bond->getBondType()) + "\n";
  }
  std::cout << "FINAL RESULT: " << std::endl;
  std::cout << bondOrders << std::endl;

  auto endTotTime = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsedTotTime = endTotTime - startTotTime;
  // std::cout << "Total time: " << elapsedTotTime.count() * 1000 << " ms"
  //           << std::endl;
  return 0;
}
