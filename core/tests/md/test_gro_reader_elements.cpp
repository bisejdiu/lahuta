#include <cstdio>
#include <memory>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "md/gro_reader.hpp"

// clang-format off
using lahuta::md::read_gro_stream_to_rwmol;

namespace {

std::string make_atom_line_no_vel(int resid, const std::string &resname, const std::string &atomname,
                                  int atomid, double x_nm, double y_nm, double z_nm) {
  char buf[256];
  std::snprintf(buf, sizeof(buf), "%5d%-5.5s%5.5s%5d%8.3f%8.3f%8.3f\n",
                resid, resname.c_str(), atomname.c_str(), atomid, x_nm, y_nm, z_nm);
  return std::string(buf);
}

} // namespace

TEST(GroReaderElements, InfersElementsAndIsotopes) {
  std::ostringstream gro;
  gro << "Elements test\n";
  gro << "7\n";
  gro << make_atom_line_no_vel(1, "NA",  "NA",  1, 0.100, 0.200, 0.300); // Ion convention: NA/NA -> Na (11)
  gro << make_atom_line_no_vel(2, "CL",  "CL",  2, 0.200, 0.300, 0.400); // Ion convention: CL/CL -> Cl (17)
  gro << make_atom_line_no_vel(3, "FE3", "FE3", 3, 0.300, 0.400, 0.500); // Curated ion: FE3/FE3 -> Fe (26)
  gro << make_atom_line_no_vel(4, "ALA", "CA",  4, 0.400, 0.500, 0.600); // Protein alpha-carbon: ALA/CA -> C (6) (must not be Ca)
  gro << make_atom_line_no_vel(5, "ALA", "HB1", 5, 0.500, 0.600, 0.700); // Protein hydrogen-like label: ALA/HB1 -> H (1)
  gro << make_atom_line_no_vel(6, "SOL", "OW",  6, 0.600, 0.700, 0.800); // Water oxygen common label: SOL/OW -> O (8)
  gro << make_atom_line_no_vel(7, "D",   "D",   7, 0.700, 0.800, 0.900); // Deuterium: D/D -> H (1) with isotope 2
  gro << "   1.00000   1.00000   1.00000\n"; // box

  std::istringstream iss(gro.str());
  auto mol = read_gro_stream_to_rwmol(iss);
  ASSERT_TRUE(mol);
  ASSERT_EQ(static_cast<int>(mol->getNumAtoms()), 7);

  // Na
  {
    const auto *a = mol->getAtomWithIdx(0);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 11);
    EXPECT_EQ(a->getIsotope(), 0);
  }

  // Cl
  {
    const auto *a = mol->getAtomWithIdx(1);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 17);
  }

  // Fe (from FE3)
  {
    const auto *a = mol->getAtomWithIdx(2);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 26);
  }

  // Protein CA should be Carbon, not Calcium
  {
    const auto *a = mol->getAtomWithIdx(3);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 6);
  }

  // HB1 -> Hydrogen
  {
    const auto *a = mol->getAtomWithIdx(4);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 1);
  }

  // OW -> Oxygen
  {
    const auto *a = mol->getAtomWithIdx(5);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 8);
  }

  // Deuterium -> Z=1, isotope=2
  {
    const auto *a = mol->getAtomWithIdx(6);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 1);
    EXPECT_EQ(a->getIsotope(), 2);
  }
}

TEST(GroReaderElements, EdgeCasesIonHeuristic) {
  std::ostringstream gro;
  gro << "Edge cases for ion heuristic\n";
  gro << "21\n";

  // Ion-like residue (curated): CLA/CL -> Cl (17)
  gro << make_atom_line_no_vel(1,  "CLA",  "CL",   1, 0.100, 0.200, 0.300);

  // Ion-like residue (curated + suffix): ZN2/ZN -> Zn (30)
  gro << make_atom_line_no_vel(2,  "ZN2",  "ZN",   2, 0.200, 0.300, 0.400);

  // Ion-like residue (curated name that isn't the pure symbol): POT/K -> K (19)
  gro << make_atom_line_no_vel(3,  "POT",  "K",    3, 0.300, 0.400, 0.500);

  // Biomolecule: ALA/ZN (two-letter element-looking atom name) -> NOT ion. Fallback: first letter 'Z' is not a valid 1-letter symbol -> unknown (0)
  gro << make_atom_line_no_vel(4,  "ALA",  "ZN",   4, 0.400, 0.500, 0.600);

  // Biomolecule: ALA/FE -> not ion-like; fallback first letter 'F' (Fluorine, 9)
  gro << make_atom_line_no_vel(5,  "ALA",  "FE",   5, 0.500, 0.600, 0.700);

  // Water-ish: SOL/OW -> O (8) via fallback first letter
  gro << make_atom_line_no_vel(6,  "SOL",  "OW",   6, 0.600, 0.700, 0.800);

  // Deuterium equality: D/D -> H (1) with isotope=2
  gro << make_atom_line_no_vel(7,  "D",    "D",    7, 0.700, 0.800, 0.900);

  // Ion-like residue via curated suffix: Ca2/CA -> Ca (20)
  gro << make_atom_line_no_vel(8,  "Ca2",  "CA",   8, 0.800, 0.900, 1.000);

  // Case + trim robustness: "  na  "/"Na" -> Na (11)
  gro << make_atom_line_no_vel(9,  "  na  ", "Na", 9, 0.900, 1.000, 1.100);

  // Ion-like via curated residue alias: NIO/NA -> Na (11)
  gro << make_atom_line_no_vel(10, "NIO",  "NA",  10, 1.000, 1.100, 1.200);

  // Single-letter deuterium in biomolecule: ALA/D -> H (1) with isotope=2 (from fallback single-letter rule)
  gro << make_atom_line_no_vel(11, "ALA",  "D",   11, 1.100, 1.200, 1.300);

  // Weird atom label with digit, non-ion residue: ALA/Z1 -> unknown (0)
  gro << make_atom_line_no_vel(12, "ALA",  "Z1",  12, 1.200, 1.300, 1.400);

  // Ion-like via trimmed equality with digits/spaces: " FE2 "/" FE3 " -> Fe (26)
  gro << make_atom_line_no_vel(13, " FE2 ",  " FE3 ",  13, 1.300, 1.400, 1.500);

  // Ion-like: zn/zn (case-insensitive) -> Zn (30)
  gro << make_atom_line_no_vel(14, "zn",     "zn",     14, 1.400, 1.500, 1.600);

  // Ion-like: "  as  "/" as " -> As (33)
  gro << make_atom_line_no_vel(15, "  as  ", " as ",   15, 1.500, 1.600, 1.700);

  // Ion-like with suffix in residue: Br2/Br -> Br (35)
  gro << make_atom_line_no_vel(16, "Br2",    "Br",     16, 1.600, 1.700, 1.800);

  // Ion-like with suffix in residue: LI3/li -> Li (3)
  gro << make_atom_line_no_vel(17, " LI3 ",  " li ",   17, 1.700, 1.800, 1.900);

  // Curated alias (3 letters): CES/CES -> Cs (55)
  gro << make_atom_line_no_vel(18, " CES ",  " CES ",  18, 1.800, 1.900, 2.000);

  // Curated alias (residue only): cal/CA -> Ca (20)
  gro << make_atom_line_no_vel(19, " cal ",  " CA ",   19, 1.900, 2.000, 2.100);

  // Water hydrogen label: SOL/HW2 -> H (1)
  gro << make_atom_line_no_vel(20, "SOL",    "HW2",    20, 2.000, 2.100, 2.200);

  // Non-ion two-letter element-looking atom under biomolecule: GLY/ZN -> unknown (0)
  gro << make_atom_line_no_vel(21, "GLY",    "ZN",     21, 2.100, 2.200, 2.300);

  // Minimal box line to finish the GRO
  gro << "   1.00000   1.00000   1.00000\n";

  std::istringstream iss(gro.str());
  auto mol = read_gro_stream_to_rwmol(iss);
  ASSERT_TRUE(mol);
  ASSERT_EQ(static_cast<int>(mol->getNumAtoms()), 21);

  // CLA/CL -> Cl
  {
    const auto *a = mol->getAtomWithIdx(0);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 17);
  }

  // ZN2/ZN -> Zn
  {
    const auto *a = mol->getAtomWithIdx(1);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 30);
  }

  // POT/K -> K
  {
    const auto *a = mol->getAtomWithIdx(2);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 19);
  }

  // ALA/ZN -> unknown (0)
  {
    const auto *a = mol->getAtomWithIdx(3);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 0);
  }

  // ALA/FE -> F (9) via first-letter fallback
  {
    const auto *a = mol->getAtomWithIdx(4);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 9);
  }

  // SOL/OW -> O (8)
  {
    const auto *a = mol->getAtomWithIdx(5);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 8);
  }

  // D/D -> H (1) + isotope 2
  {
    const auto *a = mol->getAtomWithIdx(6);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 1);
    EXPECT_EQ(a->getIsotope(), 2);
  }

  // Ca2/CA -> Ca (20)
  {
    const auto *a = mol->getAtomWithIdx(7);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 20);
  }

  // "  na  "/ "Na" -> Na  // case/trim tolerant
  {
    const auto *a = mol->getAtomWithIdx(8);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 11);
  }

  //  NIO/NA -> Na (11)
  {
    const auto *a = mol->getAtomWithIdx(9);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 11);
  }

  //  ALA/D -> H (1) + isotope 2
  {
    const auto *a = mol->getAtomWithIdx(10);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 1);
    EXPECT_EQ(a->getIsotope(), 2);
  }

  //  ALA/Z1 -> unknown (0)
  {
    const auto *a = mol->getAtomWithIdx(11);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 0);
  }

  //  Fe
  {
    const auto *a = mol->getAtomWithIdx(12);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 26);
  }

  //  Zn
  {
    const auto *a = mol->getAtomWithIdx(13);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 30);
  }

  //  As
  {
    const auto *a = mol->getAtomWithIdx(14);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 33);
  }

  //  Br
  {
    const auto *a = mol->getAtomWithIdx(15);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 35);
  }

  //  Li
  {
    const auto *a = mol->getAtomWithIdx(16);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 3);
  }

  //  Cs (from CES)
  {
    const auto *a = mol->getAtomWithIdx(17);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 55);
  }

  //  Ca (from cal alias)
  {
    const auto *a = mol->getAtomWithIdx(18);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 20);
  }

  //  H
  {
    const auto *a = mol->getAtomWithIdx(19);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 1);
  }

  // Unknown from GLY/ZN
  {
    const auto *a = mol->getAtomWithIdx(20);
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->getAtomicNum(), 0);
  }
}
