/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::ostringstream os; std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   std::copy(parts.begin(), parts.end(), std::ostream_iterator<std::string_view>(os));
 *   return os.str();
 * }();
 *
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "db/db.hpp"
#include "pipeline/ingest/lmdb.hpp"
#include "test_utils/fp_test_utils.hpp"

namespace fs = std::filesystem;
namespace lahuta::tests {
namespace P = lahuta::pipeline;

TEST(LmdbLazyLoading, PositionsViewMatchesCopyOnSampledKeys) {
  const char *env_path = std::getenv("LAHUTA_TEST_DB");
  if (!env_path) {
    GTEST_SKIP() << "LAHUTA_TEST_DB not set; skipping";
  }

  const fs::path db_path(env_path);
  if (!fs::exists(db_path)) {
    GTEST_SKIP() << "LAHUTA_TEST_DB path does not exist: " << db_path;
  }

  const std::string db_path_str(db_path.string());
  auto db = std::make_shared<LMDBDatabase>(db_path_str);

  std::vector<std::string> keys;
  db->for_each_key([&](const std::string &k) { keys.push_back(k); });
  if (keys.empty()) {
    GTEST_SKIP() << "No keys in database " << db_path_str;
  }

  const std::size_t sample_size = std::min<std::size_t>(100, keys.size());
  std::vector<std::string> sampled;
  sampled.reserve(sample_size);
  std::mt19937 rng(8351);
  std::sample(keys.begin(), keys.end(), std::back_inserter(sampled), sample_size, rng);

  for (const auto &key : sampled) {
    P::LMDBRef ref{db_path_str, /*db_name=*/"", key, db};
    auto req     = P::DataFieldSet::of({P::DataField::Positions, P::DataField::PositionsView});
    auto session = std::make_shared<P::LMDBSession>(ref, "eq_test", req);
    auto slices  = session->model_payload(req);

    ASSERT_TRUE(slices.positions) << "Positions should load positions copy for key: " << key;
    ASSERT_TRUE(slices.positions_view) << "Positions should load positions_view for key: " << key;

    const auto &copy_data = *slices.positions;
    const auto view_span  = slices.positions_view->data;
    ASSERT_EQ(copy_data.size(), view_span.size()) << "Size mismatch for key: " << key;

    for (std::size_t i = 0; i < copy_data.size(); ++i) {
      EXPECT_FLOAT_ULP_EQ(view_span[i].x, static_cast<float>(copy_data[i].x), 1, 1e-6f, i, 0);
      EXPECT_FLOAT_ULP_EQ(view_span[i].y, static_cast<float>(copy_data[i].y), 1, 1e-6f, i, 1);
      EXPECT_FLOAT_ULP_EQ(view_span[i].z, static_cast<float>(copy_data[i].z), 1, 1e-6f, i, 2);
    }
  }
}

} // namespace lahuta::tests
