#include <chrono>
#include <iostream>
#include <string>

#include "Geometry/point.h"
#include "GraphMol/RWMol.h"
#include "atom_types.hpp"
#include "contacts/halogen_bonds.hpp"
#include "contacts/hydrogen_bonds.hpp"
#include "contacts/interactions.hpp"
#include "lahuta.hpp"
#include "neighbors.hpp"
#include "nn.hpp"
#include "rings.hpp"
#include "visitor.hpp"

#include "types.hpp"

using namespace lahuta;

// TODO: 1. Remove calls to symmetrizeSSSR (DONE!, confirm before removing)

void log_ring_info(RDKit::RWMol *mol, const RingData &ring) {
  auto first_atom_index = ring.atom_ids().front();
  const auto *first_atom = mol->getAtomWithIdx(first_atom_index);
  const auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(first_atom->getMonomerInfo());
  std::cout << "Added Ring with: " << res_info->getResidueName() << " " << ring.atoms.size() << " atoms"
            << std::endl;

  auto atom_ids_copy = ring.atom_ids();
  std::sort(atom_ids_copy.begin(), atom_ids_copy.end());

  std::string atom_ids;
  for (const auto &atom_id : atom_ids_copy) {
    atom_ids += std::to_string(atom_id) + " ";
  }
  std::cout << "Atoms: " << atom_ids << std::endl;
}

void log_feature_atoms(const std::string &feature_type, const Feature &feature) {
  std::cout << "i. " << feature_type << " Charge Feature Atoms:" << std::endl;
  std::cout << "Number of atoms: " << feature.members.size() << std::endl;
  for (const auto *atom : feature.members) {
    unsigned int atom_index = atom->getIdx();
    auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    std::string atom_name = res_info->getName();
    std::string residue_name = res_info->getResidueName();
    std::string chain_id = res_info->getChainId();
    auto residue_number = res_info->getResidueNumber();

    std::cout << "i.  Atom Index: " << atom_index << ", Atom Name: " << atom_name
              << ", Residue: " << residue_name << " " << chain_id << residue_number << std::endl;
  }
}

/*bool is_cation_pi(const AtomType &f1, const AtomType &f2) {*/
/*  return (AtomTypeFlags::has(f1, AtomType::AROMATIC) && AtomTypeFlags::has(f2, AtomType::POS_IONISABLE))*/
/*         || (AtomTypeFlags::has(f1, AtomType::POS_IONISABLE) && AtomTypeFlags::has(f2, AtomType::AROMATIC));*/
/*}*/

double calculate_distance_sq(RDGeom::Point3D a, RDGeom::Point3D b) {
  double dx = a.x - b.x;
  double dy = a.y - b.y;
  double dz = a.z - b.z;
  return dx * dx + dy * dy + dz * dz;
}



int main(int argc, char const *argv[]) {
  auto start_total_time = std::chrono::high_resolution_clock::now();
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <cif file>" << std::endl;
    return 1;
  }
  std::string file_name = argv[1];

  Luni luni(file_name);
  auto mol = &luni.get_molecule();

  luni.assign_molstar_atom_types();

  const auto hydrophobic_atoms = get_atom_data(&luni, AtomType::HYDROPHOBIC);
  std::cout << "Hydrophobic Atoms: " << hydrophobic_atoms.get_data().size() << std::endl;


  Residues residues(*mol);
  EntityTypeManagerBuilder builder(*mol, residues);
  auto manager = builder.build();


  // should return const types
  AtomDataVec hp_res = manager->get_entity_type<lahuta::EntityType::Atom>(AtomType::HYDROPHOBIC);
  std::cout << "NEW Hydrophobic Atoms: " << hp_res.get_data().size() << std::endl;
  AtomDataVec hp_res_ = manager->get_entity_type<lahuta::EntityType::Atom>(AtomType::HYDROPHOBIC);
  std::cout << "NEW Hydrophobic Atoms: " << hp_res_.get_data().size() << std::endl;

  FeatureVec arom_res = manager->get_entity_type<lahuta::EntityType::Group>(AtomType::AROMATIC);
  std::cout << "AROM: " << arom_res.get_data().size() << std::endl;
  FeatureVec arom_res_ = manager->get_entity_type<lahuta::EntityType::Group>(AtomType::AROMATIC);
  std::cout << "AROM: " << arom_res_.get_data().size() << std::endl;

  std::cout << "START New contact interface" << std::endl;

  InteractionOptions opts{5.0};
  Interactions interactions(luni, opts);
  std::cout << "HBOND: \n";
  auto _1 = interactions.find_hbond_interactions();
  _1.sort_interactions();
  _1.print_interactions();
  std::cout << "WEAK HBOND: \n";
  auto _2 = interactions.find_weak_hbond_interactions();
  _2.sort_interactions();
  _2.print_interactions();
  std::cout << "HYDROPHOBIC: \n";
  auto start = std::chrono::high_resolution_clock::now();
  auto _3 = interactions.find_hydrophobic_interactions();
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << "Hydrophobic Time: " << duration.count() << " us" << std::endl;
  _3.sort_interactions();
  _3.print_interactions();
  std::cout << "HALOGEN: \n";
  auto _4 = interactions.find_halogen_interactions();
  _4.sort_interactions();
  _4.print_interactions();
  std::cout << "IONIC: \n";
  auto _5 = interactions.find_ionic_interactions();
  _5.sort_interactions();
  _5.print_interactions();
  std::cout << "METALIC: \n";
  auto _6 = interactions.find_metalic_interactions();
  _6.sort_interactions();
  _6.print_interactions();
  std::cout << "CATIONPI: \n";
  auto _7 = interactions.find_cationpi_interactions();
  _7.sort_interactions();
  _7.print_interactions();
  std::cout << "PISTACKING: \n";
  auto _8 = interactions.find_pistacking_interactions();
  _8.sort_interactions();
  _8.print_interactions();

  auto log_bond_info = [&](const RDKit::Bond *bond) {
    auto first_atom = mol->getAtomWithIdx(bond->getBeginAtomIdx());
    auto second_atom = mol->getAtomWithIdx(bond->getEndAtomIdx());
    auto *res1 = static_cast<RDKit::AtomPDBResidueInfo *>(first_atom->getMonomerInfo());
    auto *res2 = static_cast<RDKit::AtomPDBResidueInfo *>(second_atom->getMonomerInfo());
    auto atom1_name = res1->getName();
    auto atom2_name = res2->getName();
    std::string residue1 = res1->getResidueName() + "-" + std::to_string(res1->getResidueNumber());
    std::string residue2 = res2->getResidueName() + "-" + std::to_string(res2->getResidueNumber());

    auto bond_order = std::to_string(bond->getBondTypeAsDouble());
    if (bond->getIsAromatic()) {
      bond_order = "a";
    }

    std::cout << " " << residue1 << " " << residue2 << " " << bond->getBeginAtomIdx() << " "
              << bond->getEndAtomIdx() << " " << atom1_name << " " << atom2_name << " " << bond_order
              << std::endl;
  };

  int o1{}, o2{}, aromatic{};
  for (auto bond_it = mol->beginBonds(); bond_it != mol->endBonds(); ++bond_it) {
    RDKit::Bond *bond = *bond_it;

    if (bond->getBondType() == RDKit::Bond::BondType::SINGLE) {
      o1++;
    } else if (bond->getBondType() == RDKit::Bond::BondType::DOUBLE) {
      o2++;
    }
    if (bond->getIsAromatic()) {
      aromatic++;
    }
  }
  std::cout << "1: " << o1 << " 2: " << o2 << " a: " << aromatic << " t: " << o1 + o2 + aromatic << std::endl;
  std::cout << "Nr. Bonds: " << mol->getNumBonds() << std::endl;

  auto tot_time = to_ms(t() - start_total_time).count();
  std::cout << "Total Time: " << tot_time << "ms" << std::endl;

  return 0;
}
