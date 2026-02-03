/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::function<std::string(const char*, const char*, const char*)> f =
 *     [](const char* a, const char* b, const char* c) { return std::string(a) + b + c; };
 *   return f("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_INGEST_LMDB_HPP
#define LAHUTA_PIPELINE_INGEST_LMDB_HPP

#include <cstdint>
#include <cstring>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "analysis/system/records.hpp"
#include "db/model_payload.hpp"
#include "lahuta.hpp"
#include "models/topology.hpp"
#include "pipeline/data/data_requirements.hpp"
#include "pipeline/data/ingestion.hpp"
#include "pipeline/data/model_payload.hpp"
#include "pipeline/data/pipeline_item.hpp"
#include "pipeline/ingest/lmdb_coordinate_view.hpp"
#include "pipeline/ingest/lmdb_reader_slot.hpp"
#include "pipeline/session/base_session.hpp"
#include "serialization/formats.hpp"
#include "serialization/serializer.hpp"

namespace lahuta::pipeline {
namespace P = lahuta::pipeline;

class LMDBSession final : public BaseSession, public std::enable_shared_from_this<LMDBSession> {
public:
  LMDBSession(const LMDBRef &ref, std::string id, std::size_t max_inflight = 1)
      : BaseSession(std::move(id), max_inflight), key_(ref.key), db_(ref.db) {
    if (!db_) {
      throw std::invalid_argument("LMDBSession requires a valid LMDB database handle");
    }
  }

  LMDBSession(const LMDBRef &ref, std::string id, P::DataFieldSet /*unused*/, std::size_t max_inflight = 1)
      : LMDBSession(ref, std::move(id), max_inflight) {}

  const std::string &key() const noexcept { return key_; }

  P::ModelPayloadSlices model_payload(P::DataFieldSet requested) const override {
    P::ModelPayloadSlices slices;
    const bool wants_metadata       = requested.contains(P::DataField::Metadata);
    const bool wants_plddt_copy     = requested.contains(P::DataField::Plddt);
    const bool wants_plddt_view     = requested.contains(P::DataField::PlddtView);
    const bool wants_dssp_copy      = requested.contains(P::DataField::Dssp);
    const bool wants_dssp_view      = requested.contains(P::DataField::DsspView);
    const bool wants_sequence_copy  = requested.contains(P::DataField::Sequence);
    const bool wants_sequence_view  = requested.contains(P::DataField::SequenceView);
    const bool wants_frame          = requested.contains(P::DataField::Positions);
    const bool wants_positions_view = requested.contains(P::DataField::PositionsView);

    const bool needs_payload = wants_positions_view || wants_sequence_view || wants_plddt_view ||
                               wants_dssp_view || (wants_metadata && !metadata_) ||
                               (wants_plddt_copy && !plddts_) || (wants_dssp_copy && !dssp_) ||
                               (wants_sequence_copy && !sequence_);

    ReaderLease lease;
    std::optional<LoadedPayload> payload;

    if (needs_payload) lease = acquire_reader_lease(db_->get_env());

    auto ensure_payload = [&]() -> LoadedPayload & {
      if (payload) return *payload;

      if (!lease.valid()) {
        lease = acquire_reader_lease(db_->get_env());
      }

      return payload.emplace(load_payload(lease));
    };

    if (wants_metadata) {
      if (!metadata_) {
        auto &ctx = ensure_payload();
        ensure_metadata_loaded(ctx.payload);
      }
      slices.metadata = metadata_;
    }
    if (wants_plddt_copy) {
      if (!plddts_) {
        auto &ctx = ensure_payload();
        ensure_plddts_loaded(ctx.payload);
      }
      slices.plddts = plddts_;
    }
    if (wants_plddt_view) {
      auto &ctx          = ensure_payload();
      slices.plddts_view = build_plddt_view(ctx, lease);
    }
    if (wants_dssp_copy) {
      if (!dssp_) {
        auto &ctx = ensure_payload();
        ensure_dssp_loaded(ctx.payload);
      }
      slices.dssp = dssp_;
    }
    if (wants_dssp_view) {
      auto &ctx        = ensure_payload();
      slices.dssp_view = build_dssp_view(ctx, lease);
    }
    if (wants_sequence_copy) {
      if (!sequence_) {
        auto &ctx = ensure_payload();
        ensure_sequence_loaded(ctx.payload);
      }
      slices.sequence = sequence_;
    }
    if (wants_sequence_view) {
      auto &ctx            = ensure_payload();
      slices.sequence_view = build_sequence_view(ctx, lease);
    }

    std::shared_ptr<CoordinateHandle> coords_handle;
    if (wants_positions_view) {
      auto &ctx     = ensure_payload();
      coords_handle = build_coordinate_handle(ctx, lease);
    }

    if (wants_frame) {
      if (!positions_) {
        if (coords_handle) {
          positions_ = coords_handle->copy();
        } else {
          auto &ctx = ensure_payload();
          ensure_positions_loaded(ctx);
        }
      }
      slices.positions = positions_;
    }

    if (coords_handle) {
      slices.positions_view = coords_handle->view();
    }
    return slices;
  }

protected:
  std::shared_ptr<const Luni> build_system() const override {
    auto rec = ensure_record_loaded();
    if (!rec) {
      throw std::runtime_error("LMDBSession: unable to load model record for key '" + key_ + "'");
    }
    auto data_copy = rec->data; // copy to preserve original coordinates for frame consumers
    Luni luni      = Luni::from_model_data(data_copy, ModelTopologyMethod::CSR);
    return std::make_shared<Luni>(std::move(luni));
  }

private:
  struct RecordView {
    bool success = false;
    std::string_view file_path;
    std::string_view payload;
  };

  struct LoadedPayload {
    RecordView record;
    ModelPayloadView payload;

    explicit LoadedPayload(RecordView rec) : record(std::move(rec)), payload(record.payload) {}
  };

  struct CoordinateHandle {
    ReaderLease lease;
    span<const RDGeom::Point3Df> coords;
    ModelPayloadView payload;

    CoordinateHandle(ReaderLease l, span<const RDGeom::Point3Df> c, ModelPayloadView p)
        : lease(std::move(l)), coords(c), payload(std::move(p)) {}

    std::shared_ptr<const CoordinateView> view() const {
      return std::make_shared<CoordinateView>(ReaderLease(lease), coords);
    }

    std::shared_ptr<const RDGeom::POINT3D_VECT> copy() const {
      return std::make_shared<RDGeom::POINT3D_VECT>(decode_coordinates(payload));
    }
  };

  std::shared_ptr<SequenceView> build_sequence_view(const LoadedPayload &payload, ReaderLease &lease) const {
    auto seq = payload.payload.sequence();
    return std::make_shared<SequenceView>(ReaderLease(lease), seq);
  }

  std::shared_ptr<PlddtView> build_plddt_view(const LoadedPayload &payload, ReaderLease &lease) const {
    auto bytes = payload.payload.plddt_bytes();
    if (bytes.empty()) {
      return std::make_shared<PlddtView>(ReaderLease(lease), span<const pLDDTCategory>{});
    }
    if (bytes.size() % sizeof(pLDDTCategory) != 0) {
      throw std::runtime_error("LMDBSession: corrupted pLDDT payload for key '" + key_ + "'");
    }
    const auto count = bytes.size() / sizeof(pLDDTCategory);
    auto ptr         = reinterpret_cast<const pLDDTCategory *>(bytes.data());
    span<const pLDDTCategory> data(ptr, count);
    return std::make_shared<PlddtView>(ReaderLease(lease), data);
  }

  std::shared_ptr<DsspView> build_dssp_view(const LoadedPayload &payload, ReaderLease &lease) const {
    auto bytes = payload.payload.dssp_bytes();
    if (bytes.empty()) {
      return std::make_shared<DsspView>(ReaderLease(lease), span<const DSSPAssignment>{});
    }
    if (bytes.size() % sizeof(DSSPAssignment) != 0) {
      throw std::runtime_error("LMDBSession: corrupted DSSP payload for key '" + key_ + "'");
    }
    const auto count = bytes.size() / sizeof(DSSPAssignment);
    auto ptr         = reinterpret_cast<const DSSPAssignment *>(bytes.data());
    span<const DSSPAssignment> data(ptr, count);
    return std::make_shared<DsspView>(ReaderLease(lease), data);
  }

  std::shared_ptr<const analysis::ModelRecord> record() const { return ensure_record_loaded(); }

  void ensure_metadata_loaded(const ModelPayloadView &payload) const {
    std::call_once(metadata_once_, [this, &payload]() {
      auto meta     = std::make_shared<ModelMetadata>();
      auto taxonomy = payload.taxonomy_id();
      meta->ncbi_taxonomy_id.assign(taxonomy.begin(), taxonomy.end());
      auto organism = payload.organism_scientific();
      meta->organism_scientific.assign(organism.begin(), organism.end());
      metadata_ = meta;
    });
  }

  void ensure_metadata_loaded() const {
    if (metadata_) return;
    ReaderLease lease = acquire_reader_lease(db_->get_env());
    auto payload      = load_payload(lease);
    ensure_metadata_loaded(payload.payload);
  }

  void ensure_plddts_loaded(const ModelPayloadView &payload) const {
    std::call_once(plddts_once_, [this, &payload]() {
      auto bytes = payload.plddt_bytes();
      auto cats  = std::make_shared<std::vector<pLDDTCategory>>();
      if (!bytes.empty()) {
        if (bytes.size() % sizeof(pLDDTCategory) != 0) {
          throw std::runtime_error("LMDBSession: corrupted pLDDT payload for key '" + key_ + "'");
        }
        cats->resize(bytes.size() / sizeof(pLDDTCategory));
        std::memcpy(cats->data(), bytes.data(), bytes.size());
      }
      plddts_ = cats;
    });
  }

  void ensure_plddts_loaded() const {
    if (plddts_) return;
    ReaderLease lease = acquire_reader_lease(db_->get_env());
    auto payload      = load_payload(lease);
    ensure_plddts_loaded(payload.payload);
  }

  void ensure_dssp_loaded(const ModelPayloadView &payload) const {
    std::call_once(dssp_once_, [this, &payload]() {
      auto bytes       = payload.dssp_bytes();
      auto assignments = std::make_shared<std::vector<DSSPAssignment>>();
      if (!bytes.empty()) {
        if (bytes.size() % sizeof(DSSPAssignment) != 0) {
          throw std::runtime_error("LMDBSession: corrupted DSSP payload for key '" + key_ + "'");
        }
        assignments->resize(bytes.size() / sizeof(DSSPAssignment));
        std::memcpy(assignments->data(), bytes.data(), bytes.size());
      }
      dssp_ = assignments;
    });
  }

  void ensure_dssp_loaded() const {
    if (dssp_) return;
    ReaderLease lease = acquire_reader_lease(db_->get_env());
    auto payload      = load_payload(lease);
    ensure_dssp_loaded(payload.payload);
  }

  void ensure_sequence_loaded(const ModelPayloadView &payload) const {
    std::call_once(sequence_once_, [this, &payload]() {
      auto seq  = payload.sequence();
      sequence_ = std::make_shared<std::string>(seq.begin(), seq.end());
    });
  }

  void ensure_sequence_loaded() const {
    if (sequence_) return;
    ReaderLease lease = acquire_reader_lease(db_->get_env());
    auto payload      = load_payload(lease);
    ensure_sequence_loaded(payload.payload);
  }

  void ensure_positions_loaded(const LoadedPayload &payload) const {
    std::call_once(positions_once_, [this, &payload]() {
      positions_ = std::make_shared<RDGeom::POINT3D_VECT>(decode_coordinates(payload.payload));
    });
  }

  void ensure_positions_loaded() const {
    if (positions_) return;
    ReaderLease lease = acquire_reader_lease(db_->get_env());
    auto payload      = load_payload(lease);
    ensure_positions_loaded(payload);
  }

  LoadedPayload load_payload(ReaderLease &lease) const {
    if (!db_) {
      throw std::runtime_error("LMDBSession: database handle expired");
    }
    std::string_view raw;
    if (!db_->get_dbi().get(lease.txn().handle(), key_, raw)) {
      throw std::runtime_error("LMDBSession: key not found: " + key_);
    }
    auto view = parse_record_view(raw);
    if (!view.success) {
      throw std::runtime_error("LMDBSession: model record marked unsuccessful for key: " + key_);
    }
    return LoadedPayload(std::move(view));
  }

  std::shared_ptr<CoordinateHandle> build_coordinate_handle(const LoadedPayload &payload,
                                                            ReaderLease &lease) const {
    const auto count        = payload.payload.coord_count();
    const auto coords_slice = payload.payload.coord_slice();
    if (coords_slice.length == 0 && count == 0) {
      return std::make_shared<CoordinateHandle>(ReaderLease(lease),
                                                span<const RDGeom::Point3Df>{},
                                                payload.payload);
    }

    const auto payload_bytes = payload.record.payload;
    if (count == 0 || coords_slice.offset + coords_slice.length > payload_bytes.size()) {
      throw std::runtime_error("LMDBSession: corrupted coordinate payload");
    }

    const auto base                 = reinterpret_cast<const char *>(payload_bytes.data());
    const auto coords_payload_bytes = static_cast<std::size_t>(count) * CoordinateStride;
    if (coords_slice.length < coords_payload_bytes) {
      throw std::runtime_error("LMDBSession: coordinate slice too small");
    }

    //
    // Coordinate payload may not be naturally aligned because LMDB only guarantees
    // two-byte alignment for value storage. We rely on x86/ARM tolerance for
    // unaligned float access here to preserve zero-copy semantics.
    //
    auto ptr = reinterpret_cast<const RDGeom::Point3Df *>(base + coords_slice.offset);
    span<const RDGeom::Point3Df> coords_span(ptr, count);
    return std::make_shared<CoordinateHandle>(ReaderLease(lease), coords_span, payload.payload);
  }

  template <typename Fn>
  auto with_record_value(Fn &&fn, ReaderLease *lease = nullptr) const -> decltype(fn(std::string_view{})) {
    if (!db_) {
      throw std::runtime_error("LMDBSession: database handle expired");
    }
    auto local_lease       = lease ? ReaderLease{} : acquire_reader_lease(db_->get_env());
    ReaderLease &use_lease = lease ? *lease : local_lease;
    std::string_view raw;
    if (!db_->get_dbi().get(use_lease.txn().handle(), key_, raw)) {
      throw std::runtime_error("LMDBSession: key not found: " + key_);
    }
    auto result = fn(raw);
    return result;
  }

  static RecordView parse_record_view(std::string_view raw) {
    if (raw.size() < 1 + sizeof(uint32_t) * 2) {
      throw std::runtime_error("LMDBSession: corrupted record header");
    }
    RecordView view;
    const char *data   = raw.data();
    std::size_t offset = 0;
    view.success       = static_cast<unsigned char>(data[offset++]) != 0;

    uint32_t path_len = 0;
    std::memcpy(&path_len, data + offset, sizeof(path_len));
    offset += sizeof(path_len);
    if (offset + path_len + sizeof(uint32_t) > raw.size()) {
      throw std::runtime_error("LMDBSession: corrupted record path");
    }
    view.file_path = std::string_view(data + offset, path_len);
    offset += path_len;

    uint32_t blob_len = 0;
    std::memcpy(&blob_len, data + offset, sizeof(blob_len));
    offset += sizeof(blob_len);
    if (offset + blob_len > raw.size()) {
      throw std::runtime_error("LMDBSession: corrupted record payload");
    }
    view.payload = std::string_view(data + offset, blob_len);
    return view;
  }

  analysis::ModelRecord deserialize_record(std::string_view raw) const {
    return serialization::Serializer<fmt::binary, analysis::ModelRecord>::deserialize(raw.data(), raw.size());
  }

  analysis::ModelRecord load_full_record() const {
    return with_record_value(
        [this](std::string_view raw) {
          auto rec = deserialize_record(raw);
          if (!rec.success) {
            throw std::runtime_error("LMDBSession: model record marked unsuccessful for key: " + key_);
          }
          return rec;
        },
        nullptr);
  }

  std::shared_ptr<const analysis::ModelRecord> ensure_record_loaded() const {
    std::call_once(record_once_,
                   [this]() { record_ = std::make_shared<analysis::ModelRecord>(load_full_record()); });
    return record_;
  }

  static RDGeom::POINT3D_VECT decode_coordinates(const ModelPayloadView &view) {
    RDGeom::POINT3D_VECT coords;
    const auto count = view.coord_count();
    coords.resize(count);
    auto bytes = view.coords_bytes();
    if (count == 0) {
      coords.clear();
      return coords;
    }
    if (bytes.size() != count * CoordinateStride) {
      throw std::runtime_error("LMDBSession: corrupted coordinate payload");
    }
    for (std::size_t i = 0; i < count; ++i) {
      float triple[3];
      std::memcpy(triple, bytes.data() + i * sizeof(triple), sizeof(triple));
      coords[i].x = static_cast<double>(triple[0]);
      coords[i].y = static_cast<double>(triple[1]);
      coords[i].z = static_cast<double>(triple[2]);
    }
    return coords;
  }

  std::string key_;
  std::shared_ptr<LMDBDatabase> db_;

  mutable std::once_flag metadata_once_;
  mutable std::once_flag plddts_once_;
  mutable std::once_flag dssp_once_;
  mutable std::once_flag sequence_once_;
  mutable std::once_flag positions_once_;
  mutable std::once_flag record_once_;

  mutable std::shared_ptr<const ModelMetadata> metadata_;
  mutable std::shared_ptr<const std::vector<pLDDTCategory>> plddts_;
  mutable std::shared_ptr<const std::vector<DSSPAssignment>> dssp_;
  mutable std::shared_ptr<const std::string> sequence_;
  mutable std::shared_ptr<const RDGeom::POINT3D_VECT> positions_;
  mutable std::shared_ptr<analysis::ModelRecord> record_;
};

class LMDBRealizer {
public:
  void set_requirements(P::DataFieldSet req) { requirements_ = req; }

  void reset() {
    states_.clear();
    completed_.clear();
  }

  void reset(const std::string &id) {
    states_.erase(id);
    completed_.erase(id);
  }

  std::optional<PipelineItem> next(const IngestDescriptor &desc) {
    if (completed_.count(desc.id) > 0) {
      return std::nullopt;
    }

    const auto &ref = std::get<LMDBRef>(desc.origin);

    if (!needs_session()) {
      if (!completed_.insert(desc.id).second) return std::nullopt;
      PipelineItem item;
      item.session_id = desc.id;
      item.item_path  = ref.key;
      return item;
    }

    auto it = states_.find(desc.id);
    if (it == states_.end()) {
      State state;
      state.session = std::make_shared<LMDBSession>(ref, desc.id);
      it            = states_.emplace(desc.id, std::move(state)).first;
    }

    State &state = it->second;
    if (state.emitted) {
      completed_.insert(desc.id);
      states_.erase(it);
      return std::nullopt;
    }

    PipelineItem item;
    item.session_id = desc.id;
    item.item_path  = ref.key;
    item.session    = state.session;
    state.emitted   = true;
    return item;
  }

private:
  struct State {
    std::shared_ptr<LMDBSession> session;
    bool emitted = false;
  };

  bool needs_session() const { return requirements_.contains_any(P::SessionBoundFields); }

  std::unordered_map<std::string, State> states_;
  std::unordered_set<std::string> completed_;
  P::DataFieldSet requirements_ = P::DataFieldSet::none();
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_INGEST_LMDB_HPP
