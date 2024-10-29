#include <chrono>
#include <iostream>
#include <string>

#include "Geometry/point.h"
#include "GraphMol/RWMol.h"
#include "atom_types.hpp"
#include "contacts/halogen_bonds.hpp"
#include "contacts/hydrogen_bonds.hpp"
#include "contacts/interactions.hpp"
#include "contacts/utils.hpp"
#include "lahuta.hpp"
#include "neighbors.hpp"
#include "nn.hpp"
#include "nsgrid.hpp"
#include "rings.hpp"
#include "visitor.hpp"

#include "types.hpp"
/*#include "t.hpp"*/

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

bool is_same_residue(const RDKit::RWMol &mol, const RingData &ring_a, const RingData &ring_b) {
  auto atom_a = ring_a.atoms[0];
  auto atom_b = ring_b.atoms[0];
  auto info_a = static_cast<const RDKit::AtomPDBResidueInfo *>(atom_a->getMonomerInfo());
  auto info_b = static_cast<const RDKit::AtomPDBResidueInfo *>(atom_b->getMonomerInfo());
  return info_a->getResidueNumber() == info_b->getResidueNumber();
}

bool is_same_residue(const RDKit::RWMol &mol, const RDKit::Atom &atom_a, const RDKit::Atom &atom_b) {
  auto info_a = static_cast<const RDKit::AtomPDBResidueInfo *>(atom_a.getMonomerInfo());
  auto info_b = static_cast<const RDKit::AtomPDBResidueInfo *>(atom_b.getMonomerInfo());
  return info_a->getResidueNumber() == info_b->getResidueNumber();
}

bool is_cation_pi(const AtomType &f1, const AtomType &f2) {
  return (AtomTypeFlags::has(f1, AtomType::AROMATIC) && AtomTypeFlags::has(f2, AtomType::POS_IONISABLE))
         || (AtomTypeFlags::has(f1, AtomType::POS_IONISABLE) && AtomTypeFlags::has(f2, AtomType::AROMATIC));
}

double calculate_distance_sq(RDGeom::Point3D a, RDGeom::Point3D b) {
  double dx = a.x - b.x;
  double dy = a.y - b.y;
  double dz = a.z - b.z;
  return dx * dx + dy * dy + dz * dz;
}

// Project a vector onto a plane defined by a normal vector
RDGeom::Point3D project_on_plane(const RDGeom::Point3D &vector, const RDGeom::Point3D &plane_normal) {
  // subtract component along the normal
  double scalar_projection = vector.dotProduct(plane_normal);
  RDGeom::Point3D projected_vec = vector - (plane_normal * scalar_projection);
  return projected_vec;
}

double compute_in_plane_offset(
    const RDGeom::Point3D &pos_a, const RDGeom::Point3D &pos_b, const RDGeom::Point3D &normal) {
  RDGeom::Point3D vec_ab = pos_a - pos_b;
  RDGeom::Point3D projected_vec = project_on_plane(vec_ab, normal);
  double in_plane_offset = projected_vec.length();
  return in_plane_offset;
}

constexpr double angleDevMax = deg_to_rad(30.0);  // 30 degrees in radians
constexpr double deg180InRad = deg_to_rad(180.0); // π radians
constexpr double deg90InRad = deg_to_rad(90.0);   // π/2 radians
constexpr double pistacking_dist_max = 5.5;       // Maximum distance for π-stacking interactions
constexpr double offset_max = 2.1;                // Maximum offset for π-stacking interactions

void pistacking(const Luni *luni, GeometryOptions opts, Contacts &container) {

  const auto &mol = luni->get_molecule();
  const auto rings = luni->get_rings();
  /*const auto rings = luni->new_get_rings();*/

  double dist_max = 6.0;
  auto grid = FastNS(rings.centers(), dist_max);
  auto nbrs = grid.self_search();

  for (const auto &[pair, dist] : nbrs) {
    auto [ring_index_a, ring_index_b] = pair;
    const auto &ring_a = rings[ring_index_a];
    const auto &ring_b = rings[ring_index_b];

    if (is_same_residue(mol, ring_a, ring_b)) {
      continue;
    }

    auto dot_product = ring_a.norm.dotProduct(ring_b.norm);
    auto angle = std::acos(std::clamp(dot_product, -1.0, 1.0));
    if (angle > deg180InRad / 2.0) { // angle > 90 degrees
      angle = deg180InRad - angle;
    }

    double offset_a = compute_in_plane_offset(ring_a.center, ring_b.center, ring_a.norm);
    double offset_b = compute_in_plane_offset(ring_b.center, ring_a.center, ring_b.norm);

    double offset = std::min(offset_a, offset_b);

    if (offset <= offset_max) {
      if (angle <= angleDevMax) {
        std::cout << "Found PiStacking: Parallel" << std::endl;
      } else if (std::abs(angle - deg90InRad) <= angleDevMax) {
        std::cout << "Found PiStacking: T-Shaped" << std::endl;
      }
    }
  }
}

void cationpi(const Luni *luni, GeometryOptions opts, Contacts container) {

  const auto &mol = luni->get_molecule();
  const auto rings = luni->get_rings();

  auto features = get_features(luni, AtomType::POS_IONISABLE);
  if (features.features.empty()) {
    return;
  }

  double cationpi_max_dist = 6.0;
  auto grid = FastNS(rings.centers(), cationpi_max_dist);
  auto nbrs = grid.search(features.positions());

  for (const auto &[pair, _] : nbrs) {
    auto [feature_index, ring_index] = pair;
    const auto &ring = rings.rings[ring_index];
    const auto &feature = features[feature_index];

    /*auto first_ring_atom = mol.getAtomWithIdx(ring.atom_ids.front());*/
    auto first_ring_atom = ring.atoms.front();

    if (is_same_residue(mol, *first_ring_atom, *feature.members.front())) {
      continue;
    }
    auto offset = compute_in_plane_offset(feature.center, ring.center, ring.norm);

    if (offset <= 2.2) {
      std::cout << "--> Found CationPI" << std::endl;
    }
  }
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

  const auto hydrophobic_atoms = get_atom_data(&luni, AtomType::HYDROPHOBIC);
  std::cout << "Hydrophobic Atoms: " << hydrophobic_atoms.data.size() << std::endl;


  Residues residues(*mol);
  EntityTypeManagerBuilder builder(*mol, residues);
  auto manager = builder.build();

  AtomDataVec hp_res = manager->get_entity_type<lahuta::EntityType::Atom>(AtomType::HYDROPHOBIC);
  std::cout << "NEW Hydrophobic Atoms: " << hp_res.data.size() << std::endl;
  AtomDataVec hp_res_ = manager->get_entity_type<lahuta::EntityType::Atom>(AtomType::HYDROPHOBIC);
  std::cout << "NEW Hydrophobic Atoms: " << hp_res_.data.size() << std::endl;

  FeatureVec arom_res = manager->get_entity_type<lahuta::EntityType::Group>(AtomType::AROMATIC);
  std::cout << "AROM: " << arom_res.features.size() << std::endl;
  FeatureVec arom_res_ = manager->get_entity_type<lahuta::EntityType::Group>(AtomType::AROMATIC);
  std::cout << "AROM: " << arom_res_.features.size() << std::endl;

  /////////////////////////////////////////////////
  std::cout << " NEW: Pi-Stacking<  " << std::endl;
  /////////////////////////////////////////////////

  Contacts cc(&luni);
  GeometryOptions op;
  pistacking(&luni, op, cc);

  /////////////////////////////////////////////////
  std::cout << " NEW: CationPi<  " << std::endl;
  /////////////////////////////////////////////////

  Contacts cpc(&luni);
  GeometryOptions op2;
  cationpi(&luni, op2, cpc);

  std::cout << "START New contact interface" << std::endl;

  /*InteractionOptions o{5.0};*/
  /*Interactions i(&luni, o);*/
  /*auto ionic_contact = i.find_ionic_interactions();*/
  /*ionic_contact.sort_interactions();*/
  /*ionic_contact.print_interactions();*/

  // old interactions interface
  InteractionOptions opts{5.0};
  Interactions interactions(&luni, opts);
  auto _1 = interactions.find_hbond_interactions();
  auto _2 = interactions.find_weak_hbond_interactions();
  auto _3 = interactions.find_hydrophobic_interactions();
  auto _4 = interactions.find_halogen_interactions();
  auto _5 = interactions.find_ionic_interactions();
  _5.sort_interactions();
  _5.print_interactions();
  auto _6 = interactions.find_metalic_interactions();

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
