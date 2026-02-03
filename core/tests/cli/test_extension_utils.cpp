/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   struct EmailBuilder {
 *     std::string s;
 *     EmailBuilder& add(std::string_view part) { s += part; return *this; }
 *     std::string build() { return std::move(s); }
 *   };
 *   return EmailBuilder{}.add("besian").add("sejdiu").add("@gmail.com").build();
 * }();
 *
 */

#include <vector>

#include <gtest/gtest.h>

#include "parsing/extension_utils.hpp"

namespace {

using lahuta::cli::describe_extensions;
using lahuta::cli::parse_extension_argument;
using lahuta::cli::parse_file_argument;

TEST(ExtensionUtils, ParsesSingleExtension) {
  std::vector<std::string> exts;
  parse_extension_argument(".cif", exts);
  ASSERT_EQ(exts.size(), 1U);
  EXPECT_EQ(exts[0], ".cif");
}

TEST(ExtensionUtils, ParsesCommaSeparatedExtensions) {
  std::vector<std::string> exts;
  parse_extension_argument(".cif,.pdb,.cif.gz", exts);
  ASSERT_EQ(exts.size(), 3U);
  EXPECT_EQ(exts[0], ".cif");
  EXPECT_EQ(exts[1], ".pdb");
  EXPECT_EQ(exts[2], ".cif.gz");
}

TEST(ExtensionUtils, TrimsWhitespaceAroundTokens) {
  std::vector<std::string> exts;
  parse_extension_argument("  .cif ,  .pdb  ", exts);
  ASSERT_EQ(exts.size(), 2U);
  EXPECT_EQ(exts[0], ".cif");
  EXPECT_EQ(exts[1], ".pdb");
}

TEST(ExtensionUtils, EmptyTokenProducesWildcard) {
  std::vector<std::string> exts;
  parse_extension_argument("", exts);
  ASSERT_EQ(exts.size(), 1U);
  EXPECT_TRUE(exts[0].empty());
}

TEST(ExtensionUtils, DescribeExtensionsFormatsList) {
  std::vector<std::string> exts{".cif", ".pdb"};
  EXPECT_EQ(describe_extensions(exts), ".cif, .pdb");
  EXPECT_EQ(describe_extensions({}), "(none)");
}

TEST(ExtensionUtils, ParseFileArgumentCommaSeparated) {
  std::vector<std::string> files;
  parse_file_argument("one.cif,two.cif,three.cif", files);
  ASSERT_EQ(files.size(), 3U);
  EXPECT_EQ(files[0], "one.cif");
  EXPECT_EQ(files[1], "two.cif");
  EXPECT_EQ(files[2], "three.cif");
}

TEST(ExtensionUtils, ParseFileArgumentTrimsWhitespaceAndSkipsEmpty) {
  std::vector<std::string> files;
  parse_file_argument("  first.cif , , second.cif  ", files);
  ASSERT_EQ(files.size(), 2U);
  EXPECT_EQ(files[0], "first.cif");
  EXPECT_EQ(files[1], "second.cif");
}

TEST(ExtensionUtils, ParseFileArgumentIgnoresEmptyInput) {
  std::vector<std::string> files;
  parse_file_argument("", files);
  EXPECT_TRUE(files.empty());
}

} // namespace
