#ifndef LAHUTA_SOURCES_TRAJECTORY_HPP
#define LAHUTA_SOURCES_TRAJECTORY_HPP

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <functional>
#include <istream>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "lahuta.hpp"
#include "md/gro_reader.hpp"
#include "md/xtc_reader.hpp"
#include "pipeline/base_session.hpp"
#include "pipeline/frame.hpp"
#include "pipeline/ingestion.hpp"
#include "pipeline/pipeline_item.hpp"

// clang-format off
namespace {

class MemoryReadBuffer final : public std::streambuf {
public:
  MemoryReadBuffer(const char *data, std::size_t size) {
    auto *begin = const_cast<char *>(data);
    this->setg(begin, begin, begin + static_cast<std::ptrdiff_t>(size));
  }
};

} // namespace

namespace lahuta {

class TrajectorySession final : public BaseSession, public std::enable_shared_from_this<TrajectorySession> {
public:
  enum class StructureFormat { Gro, ModelCompatible };

  TrajectorySession(std::string structure, std::vector<std::string> xtcs, std::string id, std::size_t max_inflight = 128)
      : BaseSession(std::move(id), max_inflight), structure_path_(std::move(structure)),
        xtc_paths_ (std::move(xtcs)), format_(detect_format(structure_path_)) {}

  std::size_t atom_count() const {
    auto system = get_or_load_system();
    if (!system) return 0;
    return system->n_atoms();
  }

  const std::vector<std::string>& xtc_paths() const { return xtc_paths_; }

  class Frame : public FrameHandle {
  public:
    Frame(std::shared_ptr<const TrajectorySession> session, std::size_t global_idx, std::shared_ptr<const md::XtcReader::FramePacket> packet)
        : session_(std::move(session)), global_(global_idx), packet_(std::move(packet)) {}

    std::size_t index() const noexcept override { return global_; }

    std::optional<double> timestamp_ps() const noexcept override {
      return packet_ ? packet_->timestamp_ps : std::optional<double>{};
    }

    CoordinatesView load_coordinates() override {
      auto coords = get_or_materialize();
      if (!coords) {
        throw std::runtime_error("Trajectory frame coordinates unavailable");
      }
      CoordinatesView view;
      view.shared_positions = coords;
      return view;
    }

  private:
    std::shared_ptr<const RDGeom::POINT3D_VECT> get_or_materialize() const {
      std::lock_guard<std::mutex> lk(coords_mutex_);
      if (!cached_coords_) {
        cached_coords_ = materialize_locked();
      }
      return cached_coords_;
    }

    std::shared_ptr<const RDGeom::POINT3D_VECT> materialize_locked() const {
      if (!session_) {
        throw std::runtime_error("Trajectory session expired while decoding frame");
      }
      if (!packet_ || !packet_->payload) {
        throw std::runtime_error("Trajectory frame missing compressed payload");
      }

      auto buffer = session_->acquire_buffer();
      buffer->clear();
      if (packet_->atom_count > 0 && buffer->capacity() < packet_->atom_count) {
        buffer->reserve(packet_->atom_count);
      }

      MemoryReadBuffer mem_buf(packet_->payload->data(), packet_->payload->size());
      std::istream mem_stream(&mem_buf);
      md::XtcReader reader(mem_stream);
      if (!reader.read_next_into(*buffer)) {
        throw std::runtime_error("Failed to decode XTC frame from memory payload");
      }
      if (packet_->atom_count > 0 && buffer->size() != packet_->atom_count) {
        throw std::runtime_error("Decoded XTC frame atom count mismatch");
      }

      return std::shared_ptr<const RDGeom::POINT3D_VECT>(buffer, buffer.get());
    }

    std::shared_ptr<const TrajectorySession> session_;
    std::size_t global_;
    std::shared_ptr<const md::XtcReader::FramePacket> packet_;
    mutable std::mutex coords_mutex_;
    mutable std::shared_ptr<const RDGeom::POINT3D_VECT> cached_coords_;
  };

  std::shared_ptr<FrameHandle> make_frame(std::size_t global_idx,
                                          std::shared_ptr<const md::XtcReader::FramePacket> packet) const {
    if (!packet || !packet->payload) {
      throw std::invalid_argument("TrajectorySession::make_frame requires non-null payload");
    }
    return std::make_shared<Frame>(shared_from_this(), global_idx, std::move(packet));
  }

  std::unique_ptr<md::XtcReader> create_reader(std::size_t file_index) const {
    if (file_index >= xtc_paths_.size()) {
      throw std::out_of_range("TrajectorySession::create_reader index out of range");
    }
    return std::make_unique<md::XtcReader>(xtc_paths_[file_index]);
  }

protected:
  std::shared_ptr<const Luni> build_system() const override {
    switch (format_) {
      case StructureFormat::Gro: {
        auto mol = md::read_gro_to_rwmol(structure_path_);
        auto luni = Luni::create(std::move(mol));
        return std::make_shared<Luni>(std::move(luni));
      }
      case StructureFormat::ModelCompatible:
        return std::make_shared<Luni>(structure_path_);
    }
    throw std::runtime_error("TrajectorySession: unsupported structure format");
  }

private:
  friend class TrajectoryRealizer;
  static StructureFormat detect_format(const std::string& path) {
    auto to_lower  = [](unsigned char c) { return static_cast<char>(std::tolower(c)); };
    auto ends_with = [](const std::string& value, std::string_view suffix) {
      return value.size() >= suffix.size() && std::equal(suffix.rbegin(), suffix.rend(), value.rbegin());
    };

    std::string lower;
    lower.resize(path.size());
    std::transform(path.begin(), path.end(), lower.begin(), to_lower);

    return ends_with(lower, ".gro") || ends_with(lower, ".gro.gz")
               ? StructureFormat::Gro
               : StructureFormat::ModelCompatible;
  }

  std::string structure_path_;
  std::vector<std::string> xtc_paths_;
  StructureFormat format_;

  // Coordinate buffer pool, bounded by max_inflight
  std::shared_ptr<RDGeom::POINT3D_VECT> acquire_buffer() const {
    std::lock_guard<std::mutex> lk(pool_mutex_);
    if (!coord_pool_.empty()) {
      auto base = coord_pool_.back(); coord_pool_.pop_back();
      return std::shared_ptr<RDGeom::POINT3D_VECT>(
        base.get(),
        [weak_self = std::weak_ptr<const TrajectorySession>(shared_from_this()), base = std::move(base)](RDGeom::POINT3D_VECT*) mutable {
          if (auto self = weak_self.lock()) {
            self->release_buffer(std::move(base));
          }
          // else base goes out of scope and frees the buffer
        });
    }
    auto base = std::shared_ptr<RDGeom::POINT3D_VECT>(new RDGeom::POINT3D_VECT());
    return std::shared_ptr<RDGeom::POINT3D_VECT>(
      base.get(),
      [weak_self = std::weak_ptr<const TrajectorySession>(shared_from_this()), base = std::move(base)](RDGeom::POINT3D_VECT*) mutable {
        if (auto self = weak_self.lock()) {
          self->release_buffer(std::move(base));
        }
        // else base goes out of scope and frees the buffer
      });
  }

  void release_buffer(std::shared_ptr<RDGeom::POINT3D_VECT> buf) const {
    std::lock_guard<std::mutex> lk(pool_mutex_);
    if (coord_pool_.size() < max_inflight_frames()) coord_pool_.push_back(std::move(buf));
  }

  mutable std::mutex pool_mutex_;
  mutable std::vector<std::shared_ptr<RDGeom::POINT3D_VECT>> coord_pool_;
};

class TrajectoryRealizer {
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
    const auto& ref = std::get<MDRef>(desc.origin);
    if (completed_.count(desc.id) > 0) {
      return std::nullopt;
    }

    auto it = states_.find(desc.id);
    if (it == states_.end()) {
      State state;
      state.session       = std::make_shared<TrajectorySession>(ref.path, ref.xtc_paths, desc.id);
      state.file_index    = 0;
      state.global_index  = 0;
      it = states_.emplace(desc.id, std::move(state)).first;
    }

    State& state = it->second;
    auto expected_atoms = state.session->atom_count();

    while (true) {
      if (!state.reader) {
        if (state.file_index >= state.session->xtc_paths().size()) {
          completed_.insert(desc.id);
          states_.erase(it);
          return std::nullopt;
        }
        state.reader = state.session->create_reader(state.file_index);
      }

      md::XtcReader::FramePacket packet;
      if (!state.reader->read_next_packet(packet)) {
        state.reader.reset();
        ++state.file_index;
        continue;
      }

      if (!packet.payload) {
        throw std::runtime_error("XTC frame payload missing");
      }
      if (expected_atoms && packet.atom_count && packet.atom_count != expected_atoms) {
        throw std::runtime_error("XTC frame atom count does not match reference structure");
      }

      std::shared_ptr<const md::XtcReader::FramePacket> packet_ptr =
          std::make_shared<md::XtcReader::FramePacket>(std::move(packet));

      PipelineItem item;
      item.session_id   = desc.id;
      item.item_path    = ref.path;
      item.conformer_id = static_cast<std::uint64_t>(state.global_index);
      item.timestamp_ps = packet_ptr->timestamp_ps;
      item.session      = state.session;
      item.frame        = state.session->make_frame(state.global_index, std::move(packet_ptr));

      ++state.global_index;
      return item;
    }
  }

  void realize(const IngestDescriptor& desc, std::function<void(PipelineItem)> emit) {
    reset(desc.id);
    while (auto item = next(desc)) {
      emit(std::move(*item));
    }
  }

private:
  struct State {
    std::shared_ptr<TrajectorySession> session;
    std::size_t file_index   = 0;
    std::size_t global_index = 0;
    std::unique_ptr<md::XtcReader> reader;
  };

  std::unordered_map<std::string, State> states_;
  std::unordered_set<std::string> completed_;
};

} // namespace lahuta

#endif // LAHUTA_SOURCES_TRAJECTORY_HPP
