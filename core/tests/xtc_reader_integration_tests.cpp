/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto concat = [](auto&&... args) {
 *     std::string result;
 *     ((result += std::string_view(args)), ...);
 *     return result;
 *   };
 *   return concat("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

#include <array>
#include <sstream>

#include <gtest/gtest.h>

#include "md/xdr.hpp"
#include "md/xtc_reader.hpp"

using lahuta::md::XtcReader;
using lahuta::md::internal::XDRWriter;

// clang-format off
namespace {
constexpr int XtcMagic = 2023;
constexpr float NmToAngstrom = 10.0f;
} // namespace

// Write one compressed frame with byte_count=0 and a first triple that requires 3 bits.
static std::string make_truncated_xtc_frame() {
  std::ostringstream os(std::ios::binary);
  XDRWriter w(&os);

  // Header
  const int magic  = XtcMagic;
  const int natoms = 10; // > min_compressed_system_size
  const int step   = 42;
  const float time = 1.0f;
  const std::array<float, 9> box = {1.0f / NmToAngstrom, 0, 0, 0, 1.0f / NmToAngstrom, 0, 0, 0, 1.0f / NmToAngstrom};
  if (!w.write(magic))  throw std::runtime_error("Failed to write magic");
  if (!w.write(natoms)) throw std::runtime_error("Failed to write natoms");
  if (!w.write(step))   throw std::runtime_error("Failed to write step");
  if (!w.write(time))   throw std::runtime_error("Failed to write time");
  if (!w.write(box.data(), box.size())) throw std::runtime_error("Failed to write box");

  // lsize == natoms
  if (!w.write(static_cast<unsigned int>(natoms))) throw std::runtime_error("Failed to write lsize");

  // Compressed coords header
  const float precision = 0.1f;
  if (!w.write(precision)) throw std::runtime_error("Failed to write precision");

  // Bounds: minint=0, maxint=1 -> spans {2,2,2} -> bitsize = log2(8) = 3
  std::array<int, 3> minint{0, 0, 0};
  std::array<int, 3> maxint{1, 1, 1};
  if (!w.write(minint.data(), minint.size())) throw std::runtime_error("Failed to write minint");
  if (!w.write(maxint.data(), maxint.size())) throw std::runtime_error("Failed to write maxint");

  // smallidx valid (>= 9)
  const int smallidx = 9;
  if (!w.write(smallidx)) throw std::runtime_error("Failed to write smallidx");

  // byte_count (64-bit in new magic) = 0 (truncated payload)
  const std::uint64_t byte_count = 0;
  if (!w.write(byte_count)) throw std::runtime_error("Failed to write byte_count");

  return os.str();
}

TEST(XtcReaderIntegration, Compressed_TruncatedPayload_Throws) {
  const std::string frame = make_truncated_xtc_frame();
  std::istringstream is(frame, std::ios::binary);

  XtcReader reader(is);
  // The very first read should attempt to pull 3 bits for the first triple, find none, and throw.
  EXPECT_THROW({ (void)reader.read_next_frame(); }, std::exception);
}

TEST(XtcReaderIntegration, InvalidSmallIdx_Throws) {
  // Same as above but with smallidx = 8 (invalid) so header validation fails before bit IO.
  std::ostringstream os(std::ios::binary);
  XDRWriter w(&os);

  const int magic  = XtcMagic;
  const int natoms = 10;
  const int step   = 99;
  const float time = 2.0f;
  const std::array<float, 9> box = {1.0f / NmToAngstrom, 0, 0, 0, 1.0f / NmToAngstrom, 0, 0, 0, 1.0f / NmToAngstrom};
  ASSERT_TRUE(w.write(magic));
  ASSERT_TRUE(w.write(natoms));
  ASSERT_TRUE(w.write(step));
  ASSERT_TRUE(w.write(time));
  ASSERT_TRUE(w.write(box.data(), box.size()));

  ASSERT_TRUE(w.write(static_cast<unsigned int>(natoms)));
  const float precision = 0.1f;
  ASSERT_TRUE(w.write(precision));

  std::array<int, 3> minint{0, 0, 0};
  std::array<int, 3> maxint{1, 1, 1};
  ASSERT_TRUE(w.write(minint.data(), minint.size()));
  ASSERT_TRUE(w.write(maxint.data(), maxint.size()));

  const int smallidx = 8; // invalid
  ASSERT_TRUE(w.write(smallidx));

  const std::uint64_t byte_count = 0;
  ASSERT_TRUE(w.write(byte_count));

  std::istringstream is(os.str(), std::ios::binary);
  XtcReader reader(is);
  EXPECT_THROW({ (void)reader.read_next_frame(); }, std::exception);
}
