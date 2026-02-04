/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"@gmail.com", "besian", "sejdiu"};
 *   std::sort(parts.begin(), parts.end());
 *   return std::string(parts[1]) + std::string(parts[2]) + std::string(parts[0]);
 * }();
 *
 */

// Inspired by GROMACS tests (confio/trajectoryreader)

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "md/gro_reader.hpp"
#include "md/xtc_reader.hpp"

using lahuta::md::XtcReader;
using lahuta::md::internal::XDRReader;

namespace {

std::filesystem::path simdb_dir() {
  std::filesystem::path here = std::filesystem::path(__FILE__).parent_path(); // .../core/tests/md
  std::filesystem::path core = here.parent_path().parent_path();              // .../core
  return core / "data" / "simulationdatabase";
}

int gro_atom_count_from_file(const std::filesystem::path &gro) {
  std::ifstream in(gro);
  if (!in.good()) return -1;
  std::string line;
  if (!std::getline(in, line)) return -1; // title
  if (!std::getline(in, line)) return -1; // natoms
  std::istringstream iss(line);
  int n = 0;
  iss >> n;
  return n;
}

struct FirstXYZ {
  double x_ang = 0, y_ang = 0, z_ang = 0;
};

// Parse first atom xyz from GRO (nm) into Angstroms using ddist detection similar to reader
std::optional<FirstXYZ> gro_first_xyz_angstrom(const std::filesystem::path &gro) {
  std::ifstream in(gro);
  if (!in.good()) return std::nullopt;
  std::string line;
  if (!std::getline(in, line)) return std::nullopt; // title
  if (!std::getline(in, line)) return std::nullopt; // natoms
  if (!std::getline(in, line)) return std::nullopt; // first atom line
  // detect ddist by three consecutive '.'
  auto find_char = [&](std::size_t start) { return line.find('.', start); };
  auto p1 = find_char(0);
  auto p2 = (p1 != std::string::npos) ? find_char(p1 + 1) : std::string::npos;
  auto p3 = (p2 != std::string::npos) ? find_char(p2 + 1) : std::string::npos;
  if (p1 == std::string::npos || p2 == std::string::npos || p3 == std::string::npos) return std::nullopt;
  auto ddist = static_cast<int>(p2 - p1);
  if (ddist <= 0 || static_cast<int>(p3 - p2) != ddist) return std::nullopt;
  auto field = [&](std::size_t start) {
    std::size_t end = std::min(start + static_cast<std::size_t>(ddist), line.size());
    std::string s = line.substr(start, end - start);
    char *ep = nullptr;
    double v = std::strtod(s.c_str(), &ep);
    return v;
  };
  constexpr double A_per_nm = 10.0;
  const std::size_t coord_start = 20;
  FirstXYZ out;
  out.x_ang = field(coord_start + 0 * ddist) * A_per_nm;
  out.y_ang = field(coord_start + 1 * ddist) * A_per_nm;
  out.z_ang = field(coord_start + 2 * ddist) * A_per_nm;
  return out;
}

void expect_coords_finite_and_bounded(const std::vector<RDGeom::Point3D> &pts) {
  for (const auto &p : pts) {
    ASSERT_TRUE(std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z));
    // Extremely generous bound, just sanity against corrupt decoding
    ASSERT_LT(std::fabs(p.x) + std::fabs(p.y) + std::fabs(p.z), 1e9);
  }
}

struct XtcSummary {
  int frames       = 0;
  float first_time = 0.f;
  float last_time  = 0.f;
  int natoms       = 0;
};

XtcSummary scan_xtc_headers_only(const std::filesystem::path &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in.good()) {
    throw std::runtime_error("cannot open XTC: " + path.string());
  }
  XDRReader xdr(&in);
  XtcSummary out;

  while (true) {
    int magic = 0;
    if (!xdr.read(&magic)) break; // EOF cleanly
    if (magic != 1995 && magic != 2023) {
      throw std::runtime_error("XTC bad magic in scan");
    }
    int natoms = 0, step = 0;
    float time = 0.f;
    if (!xdr.read(&natoms) || !xdr.read(&step) || !xdr.read(&time)) {
      throw std::runtime_error("XTC short header");
    }
    float box[9];
    if (!xdr.read(box, 9)) {
      throw std::runtime_error("XTC short box");
    }
    if (out.frames == 0) {
      out.first_time = time;
      out.natoms = natoms;
    }
    out.last_time = time;
    ++out.frames;

    unsigned int lsize = 0;
    if (!xdr.read(&lsize)) throw std::runtime_error("XTC missing lsize");
    if (lsize <= 9) {
      // 3*lsize floats
      std::vector<float> tmp(static_cast<std::size_t>(lsize) * 3);
      if (!xdr.read(tmp.data(), tmp.size())) throw std::runtime_error("XTC short uncompressed payload");
      continue;
    }
    // compressed: precision, min/max ints, smallidx, then size and opaque bytes (with 4-byte padding)
    float precision = 0.f;
    if (!xdr.read(&precision)) throw std::runtime_error("XTC missing precision");
    int minint[3], maxint[3], smallidx;
    if (!xdr.read(minint, 3) || !xdr.read(maxint, 3) || !xdr.read(&smallidx)) {
      throw std::runtime_error("XTC missing bounds/smallidx");
    }
    std::uint64_t byte_count64 = 0;
    if (magic == 2023) {
      if (!xdr.read(reinterpret_cast<std::uint64_t *>(&byte_count64)))
        throw std::runtime_error("XTC missing size64");
    } else {
      unsigned int byte_count32 = 0;
      if (!xdr.read(&byte_count32)) throw std::runtime_error("XTC missing size32");
      byte_count64 = byte_count32;
    }
    if (byte_count64 > 0) {
      std::vector<char> skip(static_cast<std::size_t>(byte_count64));
      if (!xdr.read(skip.data(), skip.size())) throw std::runtime_error("XTC short opaque payload");
    }
  }
  return out;
}

TEST(GmxCompatGroTest, ReadsLysozymeGroBasic) {
  auto gro = simdb_dir() / "lysozyme.gro";
  int expected = gro_atom_count_from_file(gro);
  ASSERT_GT(expected, 0);
  auto mol = lahuta::md::read_gro_to_rwmol(gro.string());
  ASSERT_TRUE(mol);
  EXPECT_EQ(static_cast<int>(mol->getNumAtoms()), expected);
  auto xyz0 = gro_first_xyz_angstrom(gro);
  ASSERT_TRUE(xyz0.has_value());
  const auto &conf = mol->getConformer();
  EXPECT_NEAR(conf.getAtomPos(0).x, xyz0->x_ang, 1e-3);
  EXPECT_NEAR(conf.getAtomPos(0).y, xyz0->y_ang, 1e-3);
  EXPECT_NEAR(conf.getAtomPos(0).z, xyz0->z_ang, 1e-3);
}

struct XtcCase {
  const char *gro;
  const char *xtc;
};

std::ostream &operator<<(std::ostream &os, const XtcCase &case_) {
  return os << case_.gro << "_" << case_.xtc;
}

class GmxCompatXtcTest : public testing::TestWithParam<XtcCase> {};

TEST_P(GmxCompatXtcTest, ReadsTrajectoryFrames) {
  auto simdb = simdb_dir();
  const auto gro = simdb / GetParam().gro;
  const auto xtc = simdb / GetParam().xtc;

  int natoms_expected = gro_atom_count_from_file(gro);
  ASSERT_GT(natoms_expected, 0) << "GRO file not readable: " << gro;

  // Scan headers only to obtain exact expected reference values
  const XtcSummary hdr = scan_xtc_headers_only(xtc);
  ASSERT_EQ(hdr.natoms, natoms_expected) << "Header natoms mismatch";
  ASSERT_GE(hdr.frames, 1) << "No frames found";

  XtcReader reader(xtc.string());
  ASSERT_NO_THROW(reader.initialize());
  EXPECT_EQ(static_cast<int>(reader.natoms()), natoms_expected);
  EXPECT_TRUE(std::isfinite(reader.current_time()));
  EXPECT_FLOAT_EQ(reader.current_time(), hdr.first_time);
  expect_coords_finite_and_bounded(reader.coords());

  // Iterate remaining frames and validate
  int frames = 1;
  float prev_time = reader.current_time();
  while (reader.read_next_frame()) {
    ++frames;
    EXPECT_EQ(static_cast<int>(reader.natoms()), natoms_expected);
    EXPECT_TRUE(std::isfinite(reader.current_time()));
    EXPECT_GE(reader.current_time(), prev_time);
    prev_time = reader.current_time();
    expect_coords_finite_and_bounded(reader.coords());
  }
  EXPECT_EQ(frames, hdr.frames);
  EXPECT_FLOAT_EQ(prev_time, hdr.last_time);
}

INSTANTIATE_TEST_SUITE_P(
    SmallSystems, GmxCompatXtcTest,
    testing::Values(
        XtcCase{"spc2-traj.gro", "spc2-traj.xtc"}, XtcCase{"lysozyme.gro", "lysozyme.xtc"},
        XtcCase{"msd_coords.gro", "msd_traj.xtc"}, XtcCase{"msd_coords.gro", "msd_traj_rounding_fail.xtc"}));

// Random-access vs sequential oracle tests
struct FrameInfo {
  std::int64_t step;
  float time;
};

static std::vector<FrameInfo> build_oracle(const std::filesystem::path &xtc) {
  XtcReader r(xtc.string());
  r.initialize();
  std::vector<FrameInfo> out;
  out.push_back(FrameInfo{r.current_step(), r.current_time()});
  while (r.read_next_frame()) {
    out.push_back(FrameInfo{r.current_step(), r.current_time()});
  }
  return out;
}

static std::vector<std::size_t> sample_indices(std::size_t n, std::size_t limit = 100) {
  std::vector<std::size_t> idxs;
  if (n == 0) return idxs;
  if (n <= limit) {
    idxs.resize(n);
    std::iota(idxs.begin(), idxs.end(), 0);
    return idxs;
  }
  idxs.reserve(limit);
  for (std::size_t i = 0; i < limit; ++i) {
    std::size_t pos =
        static_cast<std::size_t>(std::llround(static_cast<long double>(i) * (n - 1) / (limit - 1)));
    idxs.push_back(pos);
  }
  return idxs;
}

TEST_P(GmxCompatXtcTest, RandomAccessByStepMatchesSequential) {
  auto simdb = simdb_dir();
  const auto xtc = simdb / GetParam().xtc;

  const auto oracle = build_oracle(xtc);
  ASSERT_FALSE(oracle.empty());

  XtcReader r(xtc.string());
  r.initialize();

  for (auto i : sample_indices(oracle.size(), 100)) {
    const auto &fi = oracle[i];
    ASSERT_TRUE(r.seek_step(fi.step)) << "seek_step failed for step=" << fi.step;
    EXPECT_EQ(r.current_step(), fi.step);
    EXPECT_FLOAT_EQ(r.current_time(), fi.time);
  }
}

TEST_P(GmxCompatXtcTest, RandomAccessByExactTimeMatchesSequential) {
  auto simdb = simdb_dir();
  const auto xtc = simdb / GetParam().xtc;

  const auto oracle = build_oracle(xtc);
  ASSERT_FALSE(oracle.empty());

  XtcReader r(xtc.string());
  r.initialize();

  for (auto i : sample_indices(oracle.size(), 100)) {
    const auto &fi = oracle[i];
    ASSERT_TRUE(r.seek_time(fi.time, /*forward_only=*/false)) << "seek_time failed for t=" << fi.time;
    EXPECT_EQ(r.current_step(), fi.step);
    EXPECT_FLOAT_EQ(r.current_time(), fi.time);
  }
}

TEST_P(GmxCompatXtcTest, RandomAccessSeekBetweenFramesLandsOnNext) {
  auto simdb = simdb_dir();
  const auto xtc = simdb / GetParam().xtc;

  const auto oracle = build_oracle(xtc);
  ASSERT_GE(oracle.size(), 2u);

  XtcReader r(xtc.string());
  r.initialize();

  // Sample up to 80 midpoints across the trajectory
  std::size_t Np = 80;
  auto idxs = sample_indices(oracle.size() - 1, Np);
  for (auto i : idxs) {
    const float t1 = oracle[i].time;
    const float t2 = oracle[i + 1].time;
    const float mid = (t1 + t2) * 0.5f;

    // Expected index: first j with time >= mid
    std::size_t exp = i + 1;
    for (std::size_t j = i; j < oracle.size(); ++j) {
      if (oracle[j].time >= mid) {
        exp = j;
        break;
      }
    }

    ASSERT_TRUE(r.seek_time(mid, /*forward_only=*/false)) << "seek_time(mid) failed at i=" << i;
    EXPECT_FLOAT_EQ(r.current_time(), oracle[exp].time);
    EXPECT_EQ(r.current_step(), oracle[exp].step);
  }
}

} // namespace
