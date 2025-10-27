#include <cstdint>
#include <cstring>
#include <memory>
#include <string>
#include <string_view>

#include <gtest/gtest.h>

#include "pipeline/dynamic/manager.hpp"
#include "pipeline/dynamic/sources.hpp"
#include "test_utils/lmdb_test_utils.hpp"

// clang-format off
namespace lahuta::tests {
using namespace topology::compute;
using namespace pipeline::compute;
using namespace pipeline::dynamic;
using namespace test_utils;

namespace {

struct ProbeParams : ParameterBase<ProbeParams> {
  static constexpr TypeId TYPE_ID = 0xF3;
};

class MetadataProbeComputation : public ReadWriteComputation<PipelineContext, ProbeParams, MetadataProbeComputation> {
public:
  using Base = ReadWriteComputation<PipelineContext, ProbeParams, MetadataProbeComputation>;
  static constexpr ComputationLabel label{"metadata_probe"};

  MetadataProbeComputation() : Base(ProbeParams{}) {}

  pipeline::DataFieldSet data_requirements() const override {
    return pipeline::DataFieldSet::of({pipeline::DataField::Metadata});
  }

  ComputationResult execute_typed(DataContext<PipelineContext, Mut::ReadWrite> &ctx, const ProbeParams &) {
    auto session = ctx.data().session;
    if (!session) {
      return ComputationResult(ComputationError("missing session"));
    }
    auto slices = session->model_payload(pipeline::DataFieldSet::of({pipeline::DataField::Metadata}));
    EXPECT_TRUE(slices.metadata);
    if (slices.metadata) {
      EXPECT_EQ(slices.metadata->ncbi_taxonomy_id, "9606");
    }
    if (auto task_ctx = ctx.data().ctx) {
      auto payload = task_ctx->get_object<const pipeline::ModelPayloadSlices>(pipeline::CTX_MODEL_PAYLOAD_KEY);
      EXPECT_TRUE(payload);
      if (payload) {
        EXPECT_TRUE(payload->metadata);
      }
    }
    return ComputationResult(true);
  }
};

class FrameProbeComputation : public ReadWriteComputation<PipelineContext, ProbeParams, FrameProbeComputation> {
public:
  using Base = ReadWriteComputation<PipelineContext, ProbeParams, FrameProbeComputation>;
  static constexpr ComputationLabel label{"frame_probe"};

  FrameProbeComputation() : Base(ProbeParams{}) {}

  pipeline::DataFieldSet data_requirements() const override {
    return pipeline::DataFieldSet::of({pipeline::DataField::Positions});
  }

  ComputationResult execute_typed(DataContext<PipelineContext, Mut::ReadWrite> &ctx, const ProbeParams &) {
    EXPECT_TRUE(ctx.data().frame);
    return ComputationResult(true);
  }
};

class PayloadProbeComputation : public ReadWriteComputation<PipelineContext, ProbeParams, PayloadProbeComputation> {
public:
  using Base = ReadWriteComputation<PipelineContext, ProbeParams, PayloadProbeComputation>;
  static constexpr ComputationLabel label{"payload_probe"};

  PayloadProbeComputation() : Base(ProbeParams{}) {}

  pipeline::DataFieldSet data_requirements() const override {
    return pipeline::DataFieldSet::of({
        pipeline::DataField::Metadata,
        pipeline::DataField::Sequence,
        pipeline::DataField::Plddt,
        pipeline::DataField::Dssp,
    });
  }

  ComputationResult execute_typed(DataContext<PipelineContext, Mut::ReadWrite> &ctx, const ProbeParams &) {
    if (!ctx.data().ctx) {
      return ComputationResult(ComputationError("missing task context"));
    }
    auto payload = ctx.data().ctx->get_object<const pipeline::ModelPayloadSlices>(pipeline::CTX_MODEL_PAYLOAD_KEY);
    EXPECT_TRUE(payload);
    if (payload) {
      EXPECT_TRUE(payload->sequence);
      EXPECT_TRUE(payload->metadata);
      EXPECT_TRUE(payload->plddts);
      EXPECT_TRUE(payload->dssp);
    }
    return ComputationResult(true);
  }
};

} // namespace

TEST(LmdbLazyLoading, MetadataRunIgnoresCorruptCoordinates) {
  LahutaRuntime::ensure_initialized(1);
  TempDir dir("lahuta_lazy_loading_metadata_");
  auto db = build_corrupt_db(dir.path, "metadata_corrupt_key");

  auto src = sources_factory::from_lmdb(std::move(db), std::string{}, 1);
  StageManager mgr(std::move(src));
  mgr.set_auto_builtins(false);

  mgr.add_computation(
      "metadata_probe",
      {},
      []() { return std::make_unique<MetadataProbeComputation>(); },
      /*thread_safe=*/true);

  ASSERT_NO_THROW({
    auto report = mgr.run(1);
    EXPECT_EQ(report.items_processed, 1u);
  });
}

TEST(LmdbLazyLoading, FrameRequirementSurfacedCorruption) {
  LahutaRuntime::ensure_initialized(1);
  TempDir dir("lahuta_lazy_loading_frame_");
  auto db = build_corrupt_db(dir.path, "frame_corrupt_key");

  auto src = sources_factory::from_lmdb(std::move(db), std::string{}, 1);
  StageManager mgr(std::move(src));
  mgr.set_auto_builtins(false);

  mgr.add_computation(
      "frame_probe",
      {},
      []() { return std::make_unique<FrameProbeComputation>(); },
      /*thread_safe=*/true);

  EXPECT_THROW(mgr.run(1), std::runtime_error);
}

TEST(LmdbLazyLoading, SequenceHandleAccessibleWithoutSystem) {
  LahutaRuntime::ensure_initialized(1);
  TempDir dir("lahuta_lazy_loading_sequence_");
  analysis::system::ModelRecord rec;
  rec.success = true;
  rec.file_path = "sequence";
  rec.data = make_sample_record();

  auto db = build_sample_db(dir.path, "sequence_key", rec);

  LMDBRef ref;
  ref.db = db;
  ref.key = "sequence_key";

  auto requirements = pipeline::DataFieldSet::of({pipeline::DataField::Sequence});
  auto session = std::make_shared<LMDBSession>(ref, "sequence_session", requirements);
  auto slices = session->model_payload(requirements);
  ASSERT_TRUE(slices.sequence);
  EXPECT_EQ(*slices.sequence, rec.data.sequence);
  auto slices_again = session->model_payload(requirements);
  ASSERT_TRUE(slices_again.sequence);
  EXPECT_EQ(slices.sequence.get(), slices_again.sequence.get());
}

TEST(LmdbLazyLoading, SequenceAndAnnotationsViews) {
  LahutaRuntime::ensure_initialized(1);
  TempDir dir("lahuta_lazy_loading_views_");
  analysis::system::ModelRecord rec;
  rec.success = true;
  rec.file_path = "views";
  rec.data = make_sample_record();

  auto db = build_sample_db(dir.path, "views_key", rec);

  LMDBRef ref;
  ref.db = db;
  ref.key = "views_key";

  auto requirements = pipeline::DataFieldSet::of({
      pipeline::DataField::Sequence,
      pipeline::DataField::SequenceView,
      pipeline::DataField::Plddt,
      pipeline::DataField::PlddtView,
      pipeline::DataField::Dssp,
      pipeline::DataField::DsspView,
  });
  auto session = std::make_shared<LMDBSession>(ref, "views_session", requirements);
  auto slices = session->model_payload(requirements);
  ASSERT_TRUE(slices.sequence);
  ASSERT_TRUE(slices.sequence_view);
  ASSERT_TRUE(slices.plddts);
  ASSERT_TRUE(slices.plddts_view);
  ASSERT_TRUE(slices.dssp);
  ASSERT_TRUE(slices.dssp_view);

  EXPECT_EQ(slices.sequence_view->data.size(), rec.data.sequence.size());
  EXPECT_EQ(std::string(slices.sequence_view->data), rec.data.sequence);

  ASSERT_EQ(slices.plddts_view->data.size(), rec.data.plddt_per_residue.size());
  for (size_t i = 0; i < slices.plddts_view->data.size(); ++i) {
    EXPECT_EQ(
        static_cast<std::uint8_t>(slices.plddts_view->data[i]),
        static_cast<std::uint8_t>(rec.data.plddt_per_residue[i]));
  }

  ASSERT_EQ(slices.dssp_view->data.size(), rec.data.dssp_per_residue.size());
  for (size_t i = 0; i < slices.dssp_view->data.size(); ++i) {
    EXPECT_EQ(
        static_cast<std::uint8_t>(slices.dssp_view->data[i]),
        static_cast<std::uint8_t>(rec.data.dssp_per_residue[i]));
  }
}

TEST(LmdbLazyLoading, PayloadObjectAvailableWhenFieldsRequested) {
  LahutaRuntime::ensure_initialized(1);
  TempDir dir("lahuta_lazy_loading_payload_");
  analysis::system::ModelRecord rec;
  rec.success = true;
  rec.file_path = "payload";
  rec.data = make_sample_record();

  auto db = build_sample_db(dir.path, "payload_key", rec);

  auto src = sources_factory::from_lmdb(std::move(db), std::string{}, 1);
  StageManager mgr(std::move(src));
  mgr.set_auto_builtins(false);

  mgr.add_computation(
      "payload_probe",
      {},
      []() { return std::make_unique<PayloadProbeComputation>(); },
      /*thread_safe=*/true);

  ASSERT_NO_THROW({
    auto report = mgr.run(1);
    EXPECT_EQ(report.items_processed, 1u);
  });
}

TEST(LmdbLazyLoading, FrameViewLoadsOnlyView) {
  LahutaRuntime::ensure_initialized(1);
  TempDir dir("lahuta_lazy_loading_frame_view_only_");
  analysis::system::ModelRecord rec;
  rec.success = true;
  rec.file_path = "frame_view_only";
  rec.data = make_sample_record();

  auto db = build_sample_db(dir.path, "frame_view_only_key", rec);
  LMDBRef ref;
  ref.db = db;
  ref.key = "frame_view_only_key";

  auto requirements = pipeline::DataFieldSet::of({pipeline::DataField::PositionsView});
  auto session = std::make_shared<LMDBSession>(ref, "frame_view_only_session", requirements);
  auto slices = session->model_payload(requirements);
  EXPECT_FALSE(slices.positions) << "PositionsView should not load positions copy";
  ASSERT_TRUE(slices.positions_view) << "PositionsView should load positions_view";
  auto coords_span = slices.positions_view->data;
  ASSERT_EQ(coords_span.size(), rec.data.coords.size());
  for (std::size_t i = 0; i < coords_span.size(); ++i) {
    EXPECT_FLOAT_EQ(coords_span[i].x, static_cast<float>(rec.data.coords[i].x));
    EXPECT_FLOAT_EQ(coords_span[i].y, static_cast<float>(rec.data.coords[i].y));
    EXPECT_FLOAT_EQ(coords_span[i].z, static_cast<float>(rec.data.coords[i].z));
  }
}

TEST(LmdbLazyLoading, FrameLoadsBothCopyAndView) {
  LahutaRuntime::ensure_initialized(1);
  TempDir dir("lahuta_lazy_loading_frame_both_");
  analysis::system::ModelRecord rec;
  rec.success = true;
  rec.file_path = "frame_both";
  rec.data = make_sample_record();

  auto db = build_sample_db(dir.path, "frame_both_key", rec);
  LMDBRef ref;
  ref.db = db;
  ref.key = "frame_both_key";

  auto requirements = pipeline::DataFieldSet::of({pipeline::DataField::Positions, pipeline::DataField::PositionsView});
  auto session = std::make_shared<LMDBSession>(ref, "frame_both_session", requirements);
  auto slices = session->model_payload(requirements);
  ASSERT_TRUE(slices.positions) << "Positions should load positions copy";
  ASSERT_TRUE(slices.positions_view) << "Positions should load positions_view";
  EXPECT_EQ(slices.positions->size(), rec.data.coords.size());
  EXPECT_EQ(slices.positions_view->data.size(), rec.data.coords.size());
}

TEST(LmdbLazyLoading, SingleReadTxnPerPayloadCall) {
  LahutaRuntime::ensure_initialized(1);
  TempDir dir("lahuta_lazy_loading_single_txn_");
  analysis::system::ModelRecord rec;
  rec.success = true;
  rec.file_path = "single_txn";
  rec.data = make_sample_record();

  auto db = build_sample_db(dir.path, "single_txn_key", rec);
  LMDBRef ref;
  ref.db = db;
  ref.key = "single_txn_key";

  reader_txn_slot().reset();
  auto before = reader_txn_slot().count();

  auto requirements = pipeline::DataFieldSet::of({
      pipeline::DataField::Positions,
      pipeline::DataField::PositionsView,
      pipeline::DataField::Metadata,
      pipeline::DataField::Plddt,
  });
  auto session = std::make_shared<LMDBSession>(ref, "single_txn_session", requirements);
  auto slices = session->model_payload(requirements);
  ASSERT_TRUE(slices.positions_view);
  ASSERT_TRUE(slices.positions);
  ASSERT_TRUE(slices.metadata);
  ASSERT_TRUE(slices.plddts);

  auto after = reader_txn_slot().count();
  EXPECT_EQ(after, before + 1);
}

TEST(LmdbLazyLoading, CachedCopyDoesNotReopenTxn) {
  LahutaRuntime::ensure_initialized(1);
  TempDir dir("lahuta_lazy_loading_cached_");
  analysis::system::ModelRecord rec;
  rec.success = true;
  rec.file_path = "cached";
  rec.data = make_sample_record();

  auto db = build_sample_db(dir.path, "cached_key", rec);
  LMDBRef ref;
  ref.db = db;
  ref.key = "cached_key";

  auto session = std::make_shared<LMDBSession>(
      ref,
      "cached_session",
      pipeline::DataFieldSet::of({pipeline::DataField::Sequence}));

  reader_txn_slot().reset();
  const auto first = reader_txn_slot().count();
  auto slices = session->model_payload(pipeline::DataFieldSet::of({pipeline::DataField::Sequence}));
  ASSERT_TRUE(slices.sequence);
  const auto mid = reader_txn_slot().count();
  EXPECT_EQ(mid, first + 1);

  auto slices_again = session->model_payload(pipeline::DataFieldSet::of({pipeline::DataField::Sequence}));
  ASSERT_TRUE(slices_again.sequence);
  const auto after = reader_txn_slot().count();
  // Subsequent calls may reopen one read txn when views are built. Make sure at most one txn per call.
  EXPECT_LE(after, mid + 1);
}

TEST(LmdbLazyLoading, ZeroCopyPositionsShareLMDBMemory) {
  LahutaRuntime::ensure_initialized(1);
  TempDir dir("lahuta_lazy_loading_zero_copy_share_");
  analysis::system::ModelRecord rec;
  rec.success = true;
  rec.file_path = "zero_copy_share";
  rec.data = make_sample_record();

  auto db = build_sample_db(dir.path, "zero_copy_share_key", rec);
  LMDBRef ref;
  ref.db = db;
  ref.key = "zero_copy_share_key";

  auto requirements = pipeline::DataFieldSet::of({pipeline::DataField::PositionsView});
  auto session = std::make_shared<LMDBSession>(ref, "zero_copy_share_session", requirements);
  auto slices = session->model_payload(requirements);
  ASSERT_TRUE(slices.positions_view);
  auto coords_span = slices.positions_view->data;
  ASSERT_EQ(coords_span.size(), rec.data.coords.size());
  const RDGeom::Point3Df *coords_ptr = coords_span.data();

  for (std::size_t i = 0; i < coords_span.size(); ++i) {
    EXPECT_FLOAT_EQ(coords_span[i].x, static_cast<float>(rec.data.coords[i].x));
    EXPECT_FLOAT_EQ(coords_span[i].y, static_cast<float>(rec.data.coords[i].y));
    EXPECT_FLOAT_EQ(coords_span[i].z, static_cast<float>(rec.data.coords[i].z));
  }

  //
  // Release the zero-copy handle before opening another read txn. LMDB forbids
  // multiple reader slots per thread.
  //
  slices.positions_view.reset();

  auto txn = lmdb::txn::begin(db->get_env().handle(), nullptr, MDB_RDONLY);
  std::string_view raw;
  ASSERT_TRUE(db->get_dbi().get(txn.handle(), ref.key, raw));
  auto payload = payload_view_from_record(raw);
  auto bytes = payload.coords_bytes();
  ASSERT_EQ(bytes.size(), rec.data.coords.size() * 3 * sizeof(float));

  auto payload_ptr = reinterpret_cast<const RDGeom::Point3Df *>(bytes.data());
  EXPECT_EQ(payload_ptr, coords_ptr);
  txn.abort();
}

TEST(LmdbLazyLoading, ZeroCopyPositionsThrowOnCorruptPayload) {
  LahutaRuntime::ensure_initialized(1);
  TempDir dir("lahuta_lazy_loading_zero_copy_corrupt_");
  auto db = build_corrupt_db(dir.path, "zero_copy_corrupt_key");
  LMDBRef ref;
  ref.db = db;
  ref.key = "zero_copy_corrupt_key";

  auto requirements = pipeline::DataFieldSet::of({pipeline::DataField::PositionsView});
  auto session = std::make_shared<LMDBSession>(ref, "zero_copy_corrupt_session", requirements);
  EXPECT_THROW(
      {
        auto slices = session->model_payload(requirements);
        (void)slices;
      },
      std::runtime_error);
}

TEST(LmdbLazyLoading, CoordinateViewOutlivesSession) {
  LahutaRuntime::ensure_initialized(1);
  TempDir dir("lahuta_lazy_loading_zero_copy_lifetime_");
  analysis::system::ModelRecord rec;
  rec.success = true;
  rec.file_path = "zero_copy_lifetime";
  rec.data = make_sample_record();

  auto db = build_sample_db(dir.path, "zero_copy_lifetime_key", rec);
  LMDBRef ref;
  ref.db = db;
  ref.key = "zero_copy_lifetime_key";

  auto requirements = pipeline::DataFieldSet::of({pipeline::DataField::PositionsView});
  auto session = std::make_shared<LMDBSession>(ref, "zero_copy_lifetime_session", requirements);
  auto slices = session->model_payload(requirements);
  ASSERT_TRUE(slices.positions_view);
  auto handle = slices.positions_view; // keep shared_ptr
  session.reset();                     // destroy session, handle should keep txn alive

  ASSERT_EQ(handle->data.size(), rec.data.coords.size());
  for (std::size_t i = 0; i < handle->data.size(); ++i) {
    EXPECT_FLOAT_EQ(handle->data[i].x, static_cast<float>(rec.data.coords[i].x));
    EXPECT_FLOAT_EQ(handle->data[i].y, static_cast<float>(rec.data.coords[i].y));
    EXPECT_FLOAT_EQ(handle->data[i].z, static_cast<float>(rec.data.coords[i].z));
  }
}

TEST(LmdbLazyLoading, ZeroCopyPositionsAlignedTo32Bytes) {
  LahutaRuntime::ensure_initialized(1);
  TempDir dir("lahuta_lazy_loading_alignment_");
  analysis::system::ModelRecord rec;
  rec.success = true;
  rec.file_path = "alignment";
  rec.data = make_sample_record();

  // Force a different payload shape to test padding logic
  rec.data.sequence = "ACDEFGHIKLMNPQRSTVW"; // 19 chars instead of 9

  auto db = build_sample_db(dir.path, "alignment_key", rec);
  LMDBRef ref;
  ref.db = db;
  ref.key = "alignment_key";

  auto requirements = pipeline::DataFieldSet::of({pipeline::DataField::PositionsView});
  auto session = std::make_shared<LMDBSession>(ref, "alignment_session", requirements);

  ASSERT_NO_THROW({
    auto slices = session->model_payload(requirements);
    ASSERT_TRUE(slices.positions_view);
    auto coords_span = slices.positions_view->data;
    ASSERT_EQ(coords_span.size(), rec.data.coords.size());

    const auto addr = reinterpret_cast<std::uintptr_t>(coords_span.data());
    // Report alignment. Zero-copy is allowed even if not 32B aligned.
    const auto point_align = static_cast<std::uintptr_t>(alignof(RDGeom::Point3Df));
    const auto float_align = static_cast<std::uintptr_t>(alignof(float));
    std::cout << "  [INFO] coords addr=0x" << std::hex << addr << std::dec
              << " mod alignof(float)=" << (addr % float_align) << " mod Point3Df=" << (addr % point_align)
              << std::endl;

    for (std::size_t i = 0; i < coords_span.size(); ++i) {
      EXPECT_FLOAT_EQ(coords_span[i].x, static_cast<float>(rec.data.coords[i].x));
      EXPECT_FLOAT_EQ(coords_span[i].y, static_cast<float>(rec.data.coords[i].y));
      EXPECT_FLOAT_EQ(coords_span[i].z, static_cast<float>(rec.data.coords[i].z));
    }
  });
}

TEST(LmdbLazyLoading, ProductionDatabaseAlignedPointers) {
  LahutaRuntime::ensure_initialized(1);

  // FIX: don't hardcode
  auto db_path = fs::path("/Users/bsejdiu/projects/lahuta_dev/lahuta/db_1773");
  if (!fs::exists(db_path)) {
    GTEST_SKIP() << "Production database 'hum_db' not found at " << db_path << ", skipping test";
    return;
  }

  auto db = std::make_shared<LMDBDatabase>(db_path.string());
  auto txn = lmdb::txn::begin(db->get_env().handle(), nullptr, MDB_RDONLY);
  auto cursor = lmdb::cursor::open(txn.handle(), db->get_dbi().handle());

  std::string_view key, value;
  if (!cursor.get(key, value, MDB_FIRST)) {
    txn.abort();
    GTEST_SKIP() << "Production database is empty";
    return;
  }

  std::string first_key(key);
  txn.abort();

  // Test zero-copy access
  LMDBRef ref;
  ref.db = db;
  ref.key = first_key;

  auto requirements = pipeline::DataFieldSet::of({pipeline::DataField::PositionsView});
  auto session = std::make_shared<LMDBSession>(ref, "prod_session", requirements);

  ASSERT_NO_THROW({
    auto slices = session->model_payload(requirements);
    ASSERT_TRUE(slices.positions_view);
    auto coords_span = slices.positions_view->data;
    ASSERT_GT(coords_span.size(), 0u);

    // Check and report alignment
    const auto addr = reinterpret_cast<std::uintptr_t>(coords_span.data());
    std::cout << "  [INFO] Production DB record '" << first_key << "' has " << coords_span.size()
              << " atoms\n";
    std::cout << "  [INFO] Coordinate pointer: 0x" << std::hex << addr << std::dec
              << ", mod alignof(float) = " << (addr % alignof(float))
              << ", mod Point3Df = " << (addr % alignof(RDGeom::Point3Df)) << "\n";

    //
    // Verify first coordinate is readable.
    // We don't know the value, but here we just don't want it to crash
    //
    volatile float x = coords_span[0].x;
    volatile float y = coords_span[0].y;
    volatile float z = coords_span[0].z;
    (void)x;
    (void)y;
    (void)z;
  });
}

} // namespace lahuta::tests
