#include <string_view>
#include <vector>

#include <gtest/gtest.h>

#include "models/parser.hpp"

namespace {

TEST(ModelParser, ParsesMarkersPastEmbeddedNulByte) {
  std::vector<char> bytes;
  auto append = [&bytes](std::string_view s) { bytes.insert(bytes.end(), s.begin(), s.end()); };

  append("data_demo\n");
  bytes.push_back('\0');
  append("_struct_ref.pdbx_seq_one_letter_code A\n");
  append("_ma_target_ref_db_details.ncbi_taxonomy_id 9606\n");
  append("_ma_target_ref_db_details.organism_scientific 'Homo sapiens'\n");
  append("ATOM ? 1.0 2.0 3.0 1.00 90.00\n");

  const auto parsed = lahuta::parse_model(bytes.data(), bytes.size());
  EXPECT_EQ(parsed.sequence, "A");
  EXPECT_EQ(parsed.metadata.ncbi_taxonomy_id, "9606");
  EXPECT_EQ(parsed.metadata.organism_scientific, "Homo sapiens");
}

} // namespace
