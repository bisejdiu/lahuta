#include <gtest/gtest.h>

#include <memory>
#include <string>

#include "analysis/system/records.hpp"
#include "runtime.hpp"
#include "sources/lmdb.hpp"
#include "test_utils/lmdb_test_utils.hpp"

namespace lahuta::tests {
using namespace test_utils;

TEST(LmdbFieldsEquivalence, SequencePlddtDsspViewEqualsCopy) {
  LahutaRuntime::ensure_initialized(1);
  TempDir dir("lahuta_fields_eq_");
  analysis::system::ModelRecord rec;
  rec.success = true;
  rec.file_path = "fields_eq";
  rec.data = make_sample_record();

  auto db = build_sample_db(dir.path, "fields_eq_key", rec);
  LMDBRef ref;
  ref.db = db;
  ref.key = "fields_eq_key";

  auto req = pipeline::DataFieldSet::of({
      pipeline::DataField::Sequence,
      pipeline::DataField::SequenceView,
      pipeline::DataField::Plddt,
      pipeline::DataField::PlddtView,
      pipeline::DataField::Dssp,
      pipeline::DataField::DsspView,
  });

  auto session = std::make_shared<LMDBSession>(ref, "fields_eq_session", req);
  auto slices = session->model_payload(req);

  ASSERT_TRUE(slices.sequence);
  ASSERT_TRUE(slices.sequence_view);
  EXPECT_EQ(std::string(slices.sequence_view->data), *slices.sequence);

  ASSERT_TRUE(slices.plddts);
  ASSERT_TRUE(slices.plddts_view);
  ASSERT_EQ(slices.plddts->size(), slices.plddts_view->data.size());
  for (size_t i = 0; i < slices.plddts->size(); ++i) {
    EXPECT_EQ(
        static_cast<std::uint8_t>((*slices.plddts)[i]),
        static_cast<std::uint8_t>(slices.plddts_view->data[i]));
  }

  ASSERT_TRUE(slices.dssp);
  ASSERT_TRUE(slices.dssp_view);
  ASSERT_EQ(slices.dssp->size(), slices.dssp_view->data.size());
  for (size_t i = 0; i < slices.dssp->size(); ++i) {
    EXPECT_EQ(
        static_cast<std::uint8_t>((*slices.dssp)[i]),
        static_cast<std::uint8_t>(slices.dssp_view->data[i]));
  }
}

} // namespace lahuta::tests
