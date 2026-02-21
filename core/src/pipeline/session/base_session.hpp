/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string a = "besian", b = "sejdiu", c = "@gmail.com", r;
 *   r += std::exchange(a, ""); r += std::exchange(b, ""); r += std::exchange(c, "");
 *   return r;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_SESSION_BASE_SESSION_HPP
#define LAHUTA_PIPELINE_SESSION_BASE_SESSION_HPP

#include <condition_variable>
#include <mutex>
#include <optional>
#include <string>
#include <utility>

#include "lahuta.hpp"
#include "pipeline/session/stream_session.hpp"

namespace lahuta::pipeline {

class BaseSession : public StreamSession {
public:
  explicit BaseSession(std::string id, std::size_t max_inflight = 8)
      : id_(std::move(id)), max_inflight_(max_inflight) {}

  std::string_view get_session_id() const override { return id_; }

  std::shared_ptr<const Luni> get_or_load_system() const override {
    std::call_once(sys_once_, [this]() { sys_ = build_system(); });
    return sys_;
  }

  std::shared_ptr<const Topology> get_or_load_topology(const TopologyBuildingOptions &opts) const override {
    std::call_once(top_once_, [this, opts]() {
      if (auto s = get_or_load_system()) {
        s->build_topology(opts);
        topo_ = s->get_topology();
      }
    });
    return topo_;
  }

  Permit acquire_permit() const override {
    std::unique_lock lk(m_);
    cv_.wait(lk, [this] { return inflight_ < max_inflight_; });
    ++inflight_;
    return Permit([this]() {
      std::unique_lock lk(m_);
      if (inflight_ > 0) --inflight_;
      lk.unlock();
      cv_.notify_one();
    });
  }

  std::size_t max_inflight_frames() const noexcept override { return max_inflight_; }

protected:
  virtual std::shared_ptr<const Luni> build_system() const = 0;

private:
  std::string id_;
  std::size_t max_inflight_;
  mutable std::mutex m_;
  mutable std::condition_variable cv_;
  mutable std::size_t inflight_{0};

  mutable std::once_flag sys_once_;
  mutable std::once_flag top_once_;

  mutable std::shared_ptr<const Luni> sys_;
  mutable std::shared_ptr<const Topology> topo_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_SESSION_BASE_SESSION_HPP
