/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::ostringstream os; std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   std::copy(parts.begin(), parts.end(), std::ostream_iterator<std::string_view>(os));
 *   return os.str();
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_INGEST_NMR_HPP
#define LAHUTA_PIPELINE_INGEST_NMR_HPP

#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <gemmi/mmread_gz.hpp>

#include "lahuta.hpp"
#include "pipeline/data/frame.hpp"
#include "pipeline/data/ingestion.hpp"
#include "pipeline/data/pipeline_item.hpp"
#include "pipeline/session/base_session.hpp"

namespace lahuta::pipeline {

class NMRSession final : public BaseSession, public std::enable_shared_from_this<NMRSession> {
public:
  explicit NMRSession(std::string path, std::string id, std::size_t max_inflight = 1)
      : BaseSession(std::move(id), max_inflight), path_(std::move(path)) {}

  struct ModelRef {
    std::size_t index;
  };

  std::size_t model_count() const {
    ensure_parsed();
    return st_.models.size();
  }

  std::vector<ModelRef> enumerate_models() const {
    ensure_parsed();
    std::vector<ModelRef> out;
    out.reserve(st_.models.size());
    for (std::size_t i = 0; i < st_.models.size(); ++i) {
      out.push_back(ModelRef{i});
    }
    return out;
  }

  class Frame : public FrameHandle {
  public:
    Frame(std::shared_ptr<const NMRSession> session, std::size_t idx)
        : session_(std::move(session)), idx_(idx) {}

    std::size_t index() const noexcept override { return idx_; }
    std::optional<double> timestamp_ps() const noexcept override { return std::nullopt; }

    CoordinatesView load_coordinates() override {
      session_->ensure_parsed();
      CoordinatesView view;
      const auto &model = session_->st_.models.at(idx_);

      auto buf = session_->acquire_buffer();
      buf->clear();
      // Estimate atoms once
      buf->reserve(session_->pooled_atom_count());

      for (const auto &chain : model.chains) {
        for (const auto &residue : chain.residues) {
          for (const auto &atom : residue.atoms) {
            buf->emplace_back(atom.pos.x, atom.pos.y, atom.pos.z);
          }
        }
      }
      view.shared_positions = std::shared_ptr<const RDGeom::POINT3D_VECT>(buf, buf.get());
      return view;
    }

  private:
    std::shared_ptr<const NMRSession> session_;
    std::size_t idx_;
  };

  std::shared_ptr<FrameHandle> make_frame(std::size_t idx) const {
    return std::make_shared<Frame>(shared_from_this(), idx);
  }

protected:
  std::shared_ptr<const Luni> build_system() const override { return std::make_shared<Luni>(path_); }

private:
  void ensure_parsed() const {
    std::call_once(parsed_once_, [this]() { st_ = gemmi::read_structure_gz(path_); });
  }

  std::string path_;
  mutable std::once_flag parsed_once_;
  mutable gemmi::Structure st_;

  // Coordinate buffer pool, bounded by max_inflight
  std::shared_ptr<RDGeom::POINT3D_VECT> acquire_buffer() const {
    std::lock_guard<std::mutex> lk(pool_mutex_);
    if (!coord_pool_.empty()) {
      auto base = coord_pool_.back();
      coord_pool_.pop_back();
      return std::shared_ptr<RDGeom::POINT3D_VECT>(
          base.get(),
          [weak_self = std::weak_ptr<const NMRSession>(shared_from_this()),
           base      = std::move(base)](RDGeom::POINT3D_VECT *) mutable {
            if (auto self = weak_self.lock()) {
              self->release_buffer(std::move(base));
            }
            // else base goes out of scope and frees the buffer
          });
    }
    // allocate new
    auto base = std::shared_ptr<RDGeom::POINT3D_VECT>(new RDGeom::POINT3D_VECT());
    return std::shared_ptr<RDGeom::POINT3D_VECT>(
        base.get(),
        [weak_self = std::weak_ptr<const NMRSession>(shared_from_this()),
         base      = std::move(base)](RDGeom::POINT3D_VECT *) mutable {
          if (auto self = weak_self.lock()) {
            self->release_buffer(std::move(base));
          }
          // else base goes out of scope and frees the buffer
        });
  }

  void release_buffer(std::shared_ptr<RDGeom::POINT3D_VECT> buf) const {
    std::lock_guard<std::mutex> lk(pool_mutex_);
    if (coord_pool_.size() < max_inflight_frames()) {
      coord_pool_.push_back(std::move(buf));
    }
  }

  std::size_t pooled_atom_count() const {
    // Try to get atom count from the first model
    std::size_t count = 0;
    if (!st_.models.empty()) {
      for (const auto &chain : st_.models.front().chains) {
        for (const auto &residue : chain.residues) {
          count += residue.atoms.size();
        }
      }
    }
    return count;
  }

  mutable std::mutex pool_mutex_;
  mutable std::vector<std::shared_ptr<RDGeom::POINT3D_VECT>> coord_pool_;
};

class NMRRealizer {
public:
  void reset() {
    states_.clear();
    completed_.clear();
  }

  void reset(const std::string &id) {
    states_.erase(id);
    completed_.erase(id);
  }

  std::optional<PipelineItem> next(const IngestDescriptor &desc) {
    const auto &ref = std::get<NMRRef>(desc.origin);
    if (completed_.count(desc.id) > 0) {
      return std::nullopt;
    }

    auto it = states_.find(desc.id);
    if (it == states_.end()) {
      auto session            = std::make_shared<NMRSession>(ref.path, desc.id);
      const auto total_models = session->model_count();
      if (total_models == 0) {
        completed_.insert(desc.id);
        return std::nullopt;
      }

      State state;
      state.session      = std::move(session);
      state.total_models = total_models;
      state.next_index   = 0;
      it                 = states_.emplace(desc.id, std::move(state)).first;
    }

    State &state = it->second;
    if (state.next_index >= state.total_models) {
      completed_.insert(desc.id);
      states_.erase(it);
      return std::nullopt;
    }

    const std::size_t model_index = state.next_index++;
    auto session                  = state.session;

    PipelineItem item;
    item.session_id   = desc.id;
    item.item_path    = ref.path;
    item.conformer_id = static_cast<std::uint64_t>(model_index);
    item.session      = session;
    item.frame        = session->make_frame(model_index);

    if (state.next_index >= state.total_models) {
      completed_.insert(desc.id);
      states_.erase(it);
    }

    return item;
  }

  void realize(const IngestDescriptor &desc, std::function<void(PipelineItem)> emit) {
    reset(desc.id);
    while (auto item = next(desc)) {
      emit(std::move(*item));
    }
  }

private:
  struct State {
    std::shared_ptr<NMRSession> session;
    std::size_t next_index   = 0;
    std::size_t total_models = 0;
  };

  std::unordered_map<std::string, State> states_;
  std::unordered_set<std::string> completed_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_INGEST_NMR_HPP
