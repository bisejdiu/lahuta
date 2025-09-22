#ifndef LAHUTA_SOURCES_LMDB_HPP
#define LAHUTA_SOURCES_LMDB_HPP

#include <memory>
#include <models/topology.hpp>
#include <mutex>
#include <optional>
#include <utility>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "analysis/system/records.hpp"
#include "db/reader.hpp"
#include "lahuta.hpp"
#include "pipeline/base_session.hpp"
#include "pipeline/ingestion.hpp"
#include "pipeline/pipeline_item.hpp"
#include "serialization/formats.hpp"

// clang-format off
namespace lahuta {

class LMDBSession final : public BaseSession, public std::enable_shared_from_this<LMDBSession> {
public:
  LMDBSession(const LMDBRef& ref, std::string id, std::size_t max_inflight = 1)
      : BaseSession(std::move(id), max_inflight), key_(ref.key), db_(ref.db) {
    if (!db_) {
      throw std::invalid_argument("LMDBSession requires a valid LMDB database handle");
    }
  }

  const std::string& key() const noexcept { return key_; }

  std::shared_ptr<FrameHandle> make_frame() const {
    auto coords = positions();
    return std::make_shared<Frame>(std::move(coords));
  }

protected:
  std::shared_ptr<const Luni> build_system() const override {
    auto rec = record();
    if (!rec) {
      throw std::runtime_error("LMDBSession: unable to load model record for key '" + key_ + "'");
    }
    auto data_copy = rec->data; // copy to preserve original coordinates for frame consumers
    Luni luni = Luni::from_model_data(data_copy, ModelTopologyMethod::CSR);
    return std::make_shared<Luni>(std::move(luni));
  }

private:
  class Frame : public FrameHandle {
  public:
    explicit Frame(std::shared_ptr<const RDGeom::POINT3D_VECT> positions)
        : positions_(std::move(positions)) {}

    std::size_t index() const noexcept override { return 0; }
    std::optional<double> timestamp_ps() const noexcept override { return std::nullopt; }

    CoordinatesView load_coordinates() override {
      CoordinatesView view;
      view.shared_positions = positions_;
      if (!view.shared_positions) {
        throw std::runtime_error("LMDBSession::Frame missing shared positions");
      }
      return view;
    }

  private:
    std::shared_ptr<const RDGeom::POINT3D_VECT> positions_;
  };

  std::shared_ptr<const analysis::system::ModelRecord> record() const {
    ensure_loaded();
    return record_;
  }

  std::shared_ptr<const RDGeom::POINT3D_VECT> positions() const {
    ensure_loaded();
    return positions_;
  }

  void ensure_loaded() const {
    std::call_once(load_once_, [this]() {
      if (!db_) {
        throw std::runtime_error("LMDBSession: database handle expired");
      }
      LMDBReader reader(db_->get_env(), db_->get_dbi());
      analysis::system::ModelRecord rec;
      if (!reader.fetch_typed<analysis::system::ModelRecord, fmt::binary>(key_, rec)) {
        throw std::runtime_error("LMDBSession: key not found: " + key_);
      }
      if (!rec.success) {
        throw std::runtime_error("LMDBSession: model record marked unsuccessful for key: " + key_);
      }
      auto positions = std::make_shared<RDGeom::POINT3D_VECT>(rec.data.coords);
      record_    = std::make_shared<analysis::system::ModelRecord>(std::move(rec));
      positions_ = std::move(positions);
    });
  }

  std::string key_;
  std::shared_ptr<LMDBDatabase> db_;

  mutable std::once_flag load_once_;
  mutable std::shared_ptr<analysis::system::ModelRecord> record_;
  mutable std::shared_ptr<const RDGeom::POINT3D_VECT> positions_;
};

class LMDBRealizer {
public:
  void reset() {
    states_   .clear();
    completed_.clear();
  }

  void reset(const std::string& id) {
    states_   .erase(id);
    completed_.erase(id);
  }

  std::optional<PipelineItem> next(const IngestDescriptor& desc) {
    if (completed_.count(desc.id) > 0) {
      return std::nullopt;
    }

    const auto& ref = std::get<LMDBRef>(desc.origin);
    auto it = states_.find(desc.id);
    if (it == states_.end()) {
      State state;
      state.session = std::make_shared<LMDBSession>(ref, desc.id);
      it = states_.emplace(desc.id, std::move(state)).first;
    }

    State& state = it->second;
    if (state.emitted) {
      completed_.insert(desc.id);
      states_.erase(it);
      return std::nullopt;
    }

    PipelineItem item;
    item.session_id = desc.id;
    item.item_path  = ref.key;
    item.session    = state.session;
    item.frame      = state.session->make_frame();
    state.emitted   = true;
    return item;
  }

private:
  struct State {
    std::shared_ptr<LMDBSession> session;
    bool emitted = false;
  };

  std::unordered_map<std::string, State> states_;
  std::unordered_set<std::string> completed_;
};

} // namespace lahuta

#endif // LAHUTA_SOURCES_LMDB_HPP
