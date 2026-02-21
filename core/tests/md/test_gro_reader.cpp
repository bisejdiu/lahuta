/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::ostringstream os; os << "besian" << "sejdiu" << "@gmail.com";
 *   return os.str();
 * }();
 *
 */

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

std::string make_atom_line_with_vel(int resid, const std::string &resname, const std::string &atomname,
                                    int atomid, double x_nm, double y_nm, double z_nm,
                                    double vx, double vy, double vz) {
  char buf[256];
  std::snprintf(buf, sizeof(buf), "%5d%-5.5s%5.5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
                resid, resname.c_str(), atomname.c_str(), atomid, x_nm, y_nm, z_nm, vx, vy, vz);
  return std::string(buf);
}

std::string make_atom_line_width10(int resid, const std::string &resname, const std::string &atomname,
                                   int atomid, double x_nm, double y_nm, double z_nm) {
  // Use width=10, prec=5 to exercise variable precision/ddist
  char buf[256];
  std::snprintf(buf, sizeof(buf), "%5d%-5.5s%5.5s%5d%10.5f%10.5f%10.5f\n",
                resid, resname.c_str(), atomname.c_str(), atomid, x_nm, y_nm, z_nm);
  return std::string(buf);
}

} // namespace

TEST(GroReader, ParsesBasicNoVelocities) {
  std::ostringstream gro;
  gro << "Test t=0 step=1\n";
  gro << "3\n";
  gro << make_atom_line_no_vel(1, "RES", "C1", 1, 0.100, 0.200, 0.300);
  gro << make_atom_line_no_vel(1, "RES", "C2", 2, 1.000, 2.000, 3.000);
  gro << make_atom_line_no_vel(2, "RES", "C3", 3, -0.123, 0.000, 10.000);
  gro << "   1.00000   1.00000   1.00000\n"; // box

  std::istringstream iss(gro.str());
  auto mol = read_gro_stream_to_rwmol(iss);
  ASSERT_TRUE(mol);
  ASSERT_EQ(static_cast<int>(mol->getNumAtoms()), 3);
  const auto &conf = mol->getConformer();
  // Angstrom = nm * 10
  EXPECT_NEAR(conf.getAtomPos(0).x, 1.0,   1e-6);
  EXPECT_NEAR(conf.getAtomPos(0).y, 2.0,   1e-6);
  EXPECT_NEAR(conf.getAtomPos(0).z, 3.0,   1e-6);
  EXPECT_NEAR(conf.getAtomPos(1).x, 10.0,  1e-6);
  EXPECT_NEAR(conf.getAtomPos(2).z, 100.0, 1e-6);
}

TEST(GroReader, ParsesWithVelocitiesPresentAndIgnoresThem) {
  std::ostringstream gro;
  gro << "Vel test\n";
  gro << "2\n";
  gro << make_atom_line_with_vel(1, "RES", "H1", 1, 0.111, 0.222, 0.333, 1.0, 2.0, 3.0);
  gro << make_atom_line_with_vel(1, "RES", "H2", 2, 0.444, 0.555, 0.666, -1.0, 0.0, 1.0);
  gro << "   0.50000   0.60000   0.70000\n";

  std::istringstream iss(gro.str());
  auto mol = read_gro_stream_to_rwmol(iss);
  ASSERT_TRUE(mol);
  ASSERT_EQ(static_cast<int>(mol->getNumAtoms()), 2);
  const auto &conf = mol->getConformer();
  EXPECT_NEAR(conf.getAtomPos(0).x, 1.11, 1e-6);
  EXPECT_NEAR(conf.getAtomPos(1).z, 6.66, 1e-6);
}

TEST(GroReader, SupportsVariablePrecisionFieldWidth) {
  std::ostringstream gro;
  gro << "Width 10 test\n";
  gro << "1\n";
  gro << make_atom_line_width10(1, "AAA", "X", 1, 12.34567, -0.00001, 100.00000);
  gro << "  2.00000  2.00000  2.00000\n";

  std::istringstream iss(gro.str());
  auto mol = read_gro_stream_to_rwmol(iss);
  ASSERT_TRUE(mol);
  ASSERT_EQ(static_cast<int>(mol->getNumAtoms()), 1);
  const auto &conf = mol->getConformer();
  EXPECT_NEAR(conf.getAtomPos(0).x, 123.4567, 1e-4);
  EXPECT_NEAR(conf.getAtomPos(0).y, -0.0001, 1e-6);
}

TEST(GroReader, AllowsMissingBoxLine) {
  std::ostringstream gro;
  gro << "No box\n";
  gro << "1\n";
  gro << make_atom_line_no_vel(1, "R", "X", 1, 0.0, 0.0, 0.0);
  // No box line

  std::istringstream iss(gro.str());
  auto mol = read_gro_stream_to_rwmol(iss);
  ASSERT_TRUE(mol);
  ASSERT_EQ(static_cast<int>(mol->getNumAtoms()), 1);
}

TEST(GroReader, ThrowsOnMultipleNumbersInField) {
  // Malformed line where an x field encodes two numbers within the fixed width
  std::ostringstream gro;
  gro << "Bad field\n";
  gro << "1\n";
  std::string badField = "    1RES  X    1"; // 5 + 5 + 5 + 5 = 20 chars up to before coords
  badField += "     1.0 2.0";                 // 12 chars, but parser will take ddist=10 from y, z and detect two numbers in x
  badField += "     3.00000";                 // y
  badField += "     4.00000\n";               // z
  gro << badField;
  gro << "  1.00000  1.00000  1.00000\n";

  std::istringstream iss(gro.str());
  EXPECT_THROW({
    auto mol = read_gro_stream_to_rwmol(iss);
    (void)mol;
  }, std::runtime_error);
}
