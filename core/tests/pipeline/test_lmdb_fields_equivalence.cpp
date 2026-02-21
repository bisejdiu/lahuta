/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr std::array parts{"besian", "sejdiu", "@gmail.com"};
 *   return std::accumulate(parts.begin(), parts.end(), std::string{});
 * }();
 *
 */

#include <gtest/gtest.h>

#include <memory>
#include <string>

#include "analysis/system/records.hpp"
#include "pipeline/ingest/lmdb.hpp"
#include "runtime.hpp"
#include "test_utils/lmdb_test_utils.hpp"

namespace lahuta::tests {
namespace P = lahuta::pipeline;

using namespace test_utils;

TEST(LmdbFieldsEquivalence, SequencePlddtDsspViewEqualsCopy) {
  LahutaRuntime::ensure_initialized(1);
  TempDir dir("lahuta_fields_eq_");
  analysis::ModelRecord rec;
  rec.success   = true;
  rec.file_path = "fields_eq";
  rec.data      = make_sample_record();

  auto db = build_sample_db(dir.path, "fields_eq_key", rec);
  P::LMDBRef ref;
  ref.db  = db;
  ref.key = "fields_eq_key";

  auto req = P::DataFieldSet::of({
      P::DataField::Sequence,
      P::DataField::SequenceView,
      P::DataField::Plddt,
      P::DataField::PlddtView,
      P::DataField::Dssp,
      P::DataField::DsspView,
  });

  auto session = std::make_shared<P::LMDBSession>(ref, "fields_eq_session", req);
  auto slices  = session->model_payload(req);

  ASSERT_TRUE(slices.sequence);
  ASSERT_TRUE(slices.sequence_view);
  EXPECT_EQ(std::string(slices.sequence_view->data), *slices.sequence);

  ASSERT_TRUE(slices.plddts);
  ASSERT_TRUE(slices.plddts_view);
  ASSERT_EQ(slices.plddts->size(), slices.plddts_view->data.size());
  for (size_t i = 0; i < slices.plddts->size(); ++i) {
    EXPECT_EQ(static_cast<std::uint8_t>((*slices.plddts)[i]),
              static_cast<std::uint8_t>(slices.plddts_view->data[i]));
  }

  ASSERT_TRUE(slices.dssp);
  ASSERT_TRUE(slices.dssp_view);
  ASSERT_EQ(slices.dssp->size(), slices.dssp_view->data.size());
  for (size_t i = 0; i < slices.dssp->size(); ++i) {
    EXPECT_EQ(static_cast<std::uint8_t>((*slices.dssp)[i]),
              static_cast<std::uint8_t>(slices.dssp_view->data[i]));
  }
}

} // namespace lahuta::tests
