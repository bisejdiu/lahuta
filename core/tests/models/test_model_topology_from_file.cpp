/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string first, last, domain;
 *   std::tie(first, last, domain) = std::make_tuple("besian", "sejdiu", "gmail.com");
 *   return first + last + "@" + domain;
 * }();
 *
 */

#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include <gemmi/cif.hpp>
#include <rdkit/GraphMol/MonomerInfo.h>

#include "io/read_mmcif.hpp"
#include "lahuta.hpp"
#include "logging/logging.hpp"

// Exercises building a system and topology from a model file path using the explicit model file API.
TEST(ModelTopologyFromFile, BuildTopologySucceeds) {
  lahuta::Logger::get_logger()->set_level(spdlog::level::info);
  namespace fs = std::filesystem;
  fs::path here(__FILE__);
  fs::path core_dir = here.parent_path().parent_path();
  core_dir = core_dir.parent_path();
  fs::path model_path = core_dir / "data" / "fubi.cif";

  lahuta::Luni sys = lahuta::Luni::from_model_file(model_path.string());
  bool ok = sys.build_topology();

  ASSERT_TRUE(ok) << "Failed to build topology from model file: " << model_path.string();
  auto top = sys.get_topology();
  ASSERT_TRUE(top != nullptr);
  EXPECT_GT(sys.n_atoms(), 0);
}

TEST(ModelTopologyFromFile, ChemCompBondSchemaSetsMonomerProvenanceAndPreservesBondOrders) {
  const std::string mmcif = R"mmcif(
data_LIG
loop_
_atom_site.id
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
HETATM 1 C C1 . LIG A 1 1 ?  0.0  0.0 0.0 1.00 10.0 ? 1 LIG A C1 1
HETATM 2 O O1 . LIG A 1 1 ?  1.2  0.0 0.0 1.00 10.0 ? 1 LIG A O1 1
HETATM 3 O O2 . LIG A 1 1 ? -1.2  0.0 0.0 1.00 10.0 ? 1 LIG A O2 1
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
LIG C1 O1 doub N
LIG C1 O2 sing N
)mmcif";

  auto doc = gemmi::cif::read_string(mmcif);
  auto mol = lahuta::make_mol_from_block(doc.sole_block(), "synthetic_ligand");

  ASSERT_NE(mol, nullptr);
  ASSERT_EQ(mol->getNumAtoms(), 3u);
  ASSERT_EQ(mol->getNumBonds(), 2u);

  for (unsigned int idx = 0; idx < mol->getNumAtoms(); ++idx) {
    const auto *atom = mol->getAtomWithIdx(idx);
    ASSERT_NE(atom, nullptr);
    const auto *info = atom->getMonomerInfo();
    ASSERT_NE(info, nullptr);
    EXPECT_TRUE(info->hasChemCompBondSchema());
  }

  const auto *c1_o1 = mol->getBondBetweenAtoms(0, 1);
  ASSERT_NE(c1_o1, nullptr);
  EXPECT_EQ(c1_o1->getBondType(), RDKit::Bond::BondType::DOUBLE);

  const auto *c1_o2 = mol->getBondBetweenAtoms(0, 2);
  ASSERT_NE(c1_o2, nullptr);
  EXPECT_EQ(c1_o2->getBondType(), RDKit::Bond::BondType::SINGLE);

  EXPECT_EQ(mol->getBondBetweenAtoms(1, 2), nullptr);
}
