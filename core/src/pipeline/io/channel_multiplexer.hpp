/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, '\0'); char* ptr = s.data();
 *   for (char c : std::string_view{"besian"}) *ptr++ = c;
 *   for (char c : std::string_view{"sejdiu"}) *ptr++ = c;
 *   *ptr++ = '@';
 *   for (char c : std::string_view{"gmail.com"}) *ptr++ = c;
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_CHANNEL_MULTIPLEXER_HPP
#define LAHUTA_PIPELINE_CHANNEL_MULTIPLEXER_HPP

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pipeline/io/backpressure.hpp"
#include "pipeline/task/emission.hpp"

namespace lahuta::pipeline {

class ChannelMultiplexer {
public:
  ChannelMultiplexer() = default;

  ChannelMultiplexer(ChannelMultiplexer &&other) noexcept {
    std::lock_guard<std::mutex> lk(other.m_);
    closed_          = other.closed_;
    ingresses_       = std::move(other.ingresses_);
    subscribers_     = std::move(other.subscribers_);
    channel_to_id_   = std::move(other.channel_to_id_);
    sink_to_ingress_ = std::move(other.sink_to_ingress_);
    next_channel_id_ = other.next_channel_id_;
  }

  ChannelMultiplexer &operator=(ChannelMultiplexer &&other) noexcept {
    if (this != &other) {
      std::scoped_lock lock(m_, other.m_);
      closed_          = other.closed_;
      ingresses_       = std::move(other.ingresses_);
      subscribers_     = std::move(other.subscribers_);
      channel_to_id_   = std::move(other.channel_to_id_);
      sink_to_ingress_ = std::move(other.sink_to_ingress_);
      next_channel_id_ = other.next_channel_id_;
    }
    return *this;
  }

  ChannelMultiplexer(const ChannelMultiplexer &)            = delete;
  ChannelMultiplexer &operator=(const ChannelMultiplexer &) = delete;

  // Connect a sink to a channel with an optional backpressure config.
  void connect(const std::string &channel, std::shared_ptr<IDynamicSink> sink,
               std::optional<BackpressureConfig> cfg = std::nullopt) {
    std::lock_guard<std::mutex> lk(m_);
    const auto ch_id = intern_channel_locked(channel);

    //
    // Defaults are applied at ingress construction time. If the same sink instance
    // is connected again later -even after defaults change - the existing ingress
    // and its captured config are reused.
    //
    BackpressureConfig effective = cfg.has_value() ? *cfg : get_default_backpressure_config();

    // one ingress per sink instance
    std::size_t idx;
    auto it = sink_to_ingress_.find(sink.get());
    if (it == sink_to_ingress_.end()) {
      idx      = ingresses_.size();
      auto ing = std::make_unique<SinkIngress>(std::move(sink), effective);
      ing->start();
      ingresses_.push_back(std::move(ing));
      sink_to_ingress_[ingresses_.back()->sink_.get()] = idx;
    } else {
      idx = it->second;
    }

    auto &subs = subscribers_[ch_id];
    if (std::find(subs.begin(), subs.end(), idx) == subs.end()) subs.push_back(idx);
  }

  // Producer-side API: capture payload by move and enqueue views to subscribed sinks
  void emit(Emission &&e) {
    std::vector<std::size_t> subs;
    uint32_t ch_id = 0;
    {
      std::lock_guard<std::mutex> lk(m_);
      if (closed_) return; // ignore after close
      ch_id   = intern_channel_locked(e.channel);
      auto it = subscribers_.find(ch_id);
      if (it != subscribers_.end()) subs = it->second;
    }
    if (subs.empty()) return;

    if (subs.size() == 1) { // optimization for single-subscriber case
      auto buf      = std::make_shared<std::string>(std::move(e.payload));
      const auto sz = buf->size();
      auto &ing     = *ingresses_[subs[0]];
      (void)ing.offer(ch_id, std::move(buf), sz);
      return;
    }

    auto buf      = std::make_shared<std::string>(std::move(e.payload));
    const auto sz = buf->size();
    for (auto idx : subs) {
      auto &ing = *ingresses_[idx];
      (void)ing.offer(ch_id, buf, sz);
    }
  }

  // Stop ingress and wait for writers with a deadline. Throws on required sink failure/timeouts.
  void close_and_flush(std::chrono::milliseconds timeout) {
    std::vector<SinkIngress *> ing;
    {
      std::lock_guard<std::mutex> lk(m_);
      if (closed_) return;
      closed_ = true;
      for (auto &p : ingresses_) {
        p->close_queue();
        ing.push_back(p.get());
      }
    }

    const auto deadline      = std::chrono::steady_clock::now() + timeout;
    bool any_required_failed = false;
    std::string err;
    for (auto *si : ing) {
      const bool ok = si->join_until(deadline);
      if (!ok) {
        // couldn't finish in time
        if (si->cfg_.required) {
          any_required_failed = true;
          err                 = "sink deadline exceeded";
        }
      } else if (si->failed_ && si->cfg_.required) {
        any_required_failed = true;
        err                 = si->error_;
      }
    }
    if (any_required_failed) {
      throw std::runtime_error("ChannelMultiplexer.close_and_flush: required sink failed: " + err);
    }
  }

  // Prepare for a new run: if previously closed, restart all writer threads and reopen ingress queues.
  void reopen_if_closed() {
    std::lock_guard<std::mutex> lk(m_);
    if (!closed_) return;
    closed_ = false;
    for (auto &p : ingresses_) {
      // Recreate ingress with the same sink and config. Start a fresh writer thread
      auto new_ing = std::make_unique<SinkIngress>(p->sink_, p->cfg_);
      new_ing->start();
      p.swap(new_ing);
    }
  }

  struct SinkStatsSnapshot {
    std::string sink_name;
    uint64_t enqueued_msgs  = 0;
    uint64_t enqueued_bytes = 0;
    uint64_t written_msgs   = 0;
    uint64_t written_bytes  = 0;
    uint64_t stalled_ns     = 0;
    uint64_t drops          = 0;
    std::size_t queue_msgs  = 0;
    std::size_t queue_bytes = 0;

    std::size_t writer_threads         = 0;
    std::size_t queue_high_water_msgs  = 0;
    std::size_t queue_high_water_bytes = 0;
    std::size_t active_writers         = 0;
  };

  std::vector<SinkStatsSnapshot> stats() const {
    std::vector<SinkStatsSnapshot> out;
    std::lock_guard<std::mutex> lk(m_);
    out.reserve(ingresses_.size());
    for (auto const &p : ingresses_) {
      const auto *si = p.get();
      SinkStatsSnapshot s;
      s.sink_name              = {};
      s.enqueued_msgs          = si->enq_msgs_.load(std::memory_order_relaxed);
      s.enqueued_bytes         = si->enq_bytes_.load(std::memory_order_relaxed);
      s.written_msgs           = si->written_msgs_.load(std::memory_order_relaxed);
      s.written_bytes          = si->written_bytes_.load(std::memory_order_relaxed);
      s.stalled_ns             = si->stall_ns_accum_.load(std::memory_order_relaxed);
      s.drops                  = si->drops_counter_.load(std::memory_order_relaxed);
      s.queue_msgs             = si->queue_.size_msgs();
      s.queue_bytes            = si->queue_.size_bytes();
      s.writer_threads         = si->cfg_.writer_threads;
      s.queue_high_water_msgs  = si->queue_.high_water_msgs();
      s.queue_high_water_bytes = si->queue_.high_water_bytes();
      s.active_writers         = si->active_writers_.load(std::memory_order_relaxed);
      out.push_back(std::move(s));
    }
    return out;
  }

private:
  uint32_t intern_channel_locked(const std::string &ch) {
    auto it = channel_to_id_.find(ch);
    if (it != channel_to_id_.end()) return it->second;
    uint32_t id        = next_channel_id_++;
    channel_to_id_[ch] = id;
    return id;
  }

private:
  mutable std::mutex m_;
  bool closed_ = false;

  std::vector<std::unique_ptr<SinkIngress>> ingresses_;
  std::unordered_map<IDynamicSink *, std::size_t> sink_to_ingress_;
  std::unordered_map<uint32_t, std::vector<std::size_t>> subscribers_;

  // channel interning
  std::unordered_map<std::string, uint32_t> channel_to_id_;
  uint32_t next_channel_id_ = 0;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_CHANNEL_MULTIPLEXER_HPP
