/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"@gmail.com", "besian", "sejdiu"};
 *   std::sort(parts.begin(), parts.end());
 *   return std::string(parts[1]) + std::string(parts[2]) + std::string(parts[0]);
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_BACKPRESSURE_HPP
#define LAHUTA_PIPELINE_BACKPRESSURE_HPP

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

#include "pipeline/io/sink_iface.hpp"
#include "pipeline/task/emission.hpp"

namespace lahuta::pipeline {

// Backpressure policy when a per-sink queue is full
enum class OnFull { Block, DropLatest, DropOldest };

struct BackpressureConfig {
  std::size_t max_queue_msgs  = 65536;
  std::size_t max_queue_bytes = 1ull * 1024 * 1024 * 1024; // = 1 GiB
  std::size_t max_batch_msgs  = 256;
  std::size_t max_batch_bytes = 50 * 1024 * 1024; // = 50 MiB
  std::size_t writer_threads  = 1;
  // Producer wait-slice while blocking on full queue. Not a hard timeout.
  // Total blocking under OnFull::Block is not bounded. This is the maximum duration
  // of a single wait before re-checking capacity. With correct (?) condition-variable
  // notifications, producers wake immediately when space is freed.
  // Low-latency 5-20 ms, balanced default 100 ms, background 200-500 ms.
  std::chrono::milliseconds offer_wait_slice{100};
  OnFull on_full = OnFull::Block;
  bool required  = true; // failure aborts pipeline
};

inline void validate_config(const BackpressureConfig &cfg) {
  if (cfg.max_queue_msgs == 0) throw std::invalid_argument("max_queue_msgs  must be > 0");
  if (cfg.max_queue_bytes == 0) throw std::invalid_argument("max_queue_bytes must be > 0");
  if (cfg.max_batch_msgs == 0) throw std::invalid_argument("max_batch_msgs  must be > 0");
  if (cfg.max_batch_bytes == 0) throw std::invalid_argument("max_batch_bytes must be > 0");
  if (cfg.writer_threads == 0) throw std::invalid_argument("writer_threads  must be > 0");
  if (cfg.max_batch_bytes > cfg.max_queue_bytes) {
    throw std::invalid_argument("max_batch_bytes cannot exceed max_queue_bytes");
  }
}

namespace detail {
inline std::mutex &default_cfg_mutex() {
  static std::mutex m;
  return m;
}

inline BackpressureConfig &default_cfg_ref() {
  static BackpressureConfig cfg;
  return cfg;
}
} // namespace detail

inline BackpressureConfig get_default_backpressure_config() {
  std::lock_guard<std::mutex> lk(detail::default_cfg_mutex());
  return detail::default_cfg_ref();
}

inline void set_default_backpressure_config(const BackpressureConfig &cfg) {
  validate_config(cfg);
  std::lock_guard<std::mutex> lk(detail::default_cfg_mutex());
  detail::default_cfg_ref() = cfg;
}

inline void set_default_max_queue_bytes(std::size_t bytes) {
  auto cfg            = get_default_backpressure_config();
  cfg.max_queue_bytes = bytes;
  if (cfg.max_batch_bytes > bytes) cfg.max_batch_bytes = bytes;
  set_default_backpressure_config(cfg);
}

// Node stored in per-sink queue
struct QueueNode {
  uint32_t channel_id = 0;
  std::shared_ptr<const std::string> payload;
  std::size_t size = 0;
};

class BoundedQueue {
public:
  BoundedQueue(std::size_t max_msgs, std::size_t max_bytes) : max_msgs_(max_msgs), max_bytes_(max_bytes) {}

  bool offer(QueueNode n, const BackpressureConfig &cfg, std::atomic<uint64_t> *stall_ns_accum,
             std::atomic<uint64_t> *drops_counter) {
    const auto need = n.size;
    if (need > max_bytes_) {
      if (drops_counter) ++(*drops_counter);
      return false;
    }
    std::unique_lock<std::mutex> lk(m_);
    for (;;) {
      if (closed_) return false;
      if ((q_.size() < max_msgs_) && (bytes_ + need <= max_bytes_)) {
        q_.push_back(std::move(n));
        bytes_ += need;
        if (q_.size() > high_water_msgs_) high_water_msgs_ = q_.size();
        if (bytes_ > high_water_bytes_) high_water_bytes_ = bytes_;
        lk.unlock();
        cv_not_empty_.notify_one();
        return true;
      }
      switch (cfg.on_full) {
        case OnFull::DropLatest:
          if (drops_counter) ++(*drops_counter);
          return false;
        case OnFull::DropOldest:
          if (!q_.empty()) {
            bytes_ -= q_.front().size;
            q_.pop_front();
            if (drops_counter) ++(*drops_counter);
            continue;
          }
          [[fallthrough]];
        case OnFull::Block: {
          auto t0 = std::chrono::steady_clock::now();
          cv_not_full_.wait_for(lk, cfg.offer_wait_slice, [&] {
            return closed_ || ((q_.size() < max_msgs_) && (bytes_ + need <= max_bytes_));
          });
          if (closed_) return false;
          auto t1 = std::chrono::steady_clock::now();
          if (stall_ns_accum)
            *stall_ns_accum += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
          continue;
        }
      }
    }
  }

  bool pop_batch(std::vector<QueueNode> &out, std::size_t max_msgs, std::size_t max_bytes) {
    std::unique_lock<std::mutex> lk(m_);
    cv_not_empty_.wait(lk, [&] { return closed_ || !q_.empty(); });
    if (q_.empty() && closed_) return false;
    std::size_t taken_msgs = 0, taken_bytes = 0;
    while (!q_.empty() && taken_msgs < max_msgs) {
      const auto front_size = q_.front().size;
      // Guarantee forward progress: always take the first item even if it exceeds max_bytes.
      if (taken_msgs > 0 && taken_bytes + front_size > max_bytes) break;
      QueueNode node = std::move(q_.front());
      q_.pop_front();
      bytes_ -= node.size;
      taken_bytes += node.size;
      ++taken_msgs;
      out.push_back(std::move(node));
    }
    const auto freed_msgs = taken_msgs;
    lk.unlock();
    // Selective wake-up: notify_one for a single freed slot, notify_all for multiple.
    if (freed_msgs > 1) {
      cv_not_full_.notify_all();
    } else if (freed_msgs == 1) {
      cv_not_full_.notify_one();
    }
    return !out.empty();
  }

  void close() {
    std::lock_guard<std::mutex> lk(m_);
    closed_ = true;
    cv_not_empty_.notify_all();
    cv_not_full_.notify_all();
  }

  std::size_t size_msgs() const {
    std::lock_guard<std::mutex> lk(m_);
    return q_.size();
  }

  std::size_t size_bytes() const {
    std::lock_guard<std::mutex> lk(m_);
    return bytes_;
  }

  std::size_t high_water_msgs() const {
    std::lock_guard<std::mutex> lk(m_);
    return high_water_msgs_;
  }

  std::size_t high_water_bytes() const {
    std::lock_guard<std::mutex> lk(m_);
    return high_water_bytes_;
  }

private:
  std::size_t max_msgs_  = 0;
  std::size_t max_bytes_ = 0;
  std::deque<QueueNode> q_;
  std::size_t bytes_            = 0;
  std::size_t high_water_msgs_  = 0;
  std::size_t high_water_bytes_ = 0;
  mutable std::mutex m_;
  std::condition_variable cv_not_full_;
  std::condition_variable cv_not_empty_;
  bool closed_ = false;
};

struct SinkIngress {
  explicit SinkIngress(std::shared_ptr<IDynamicSink> sink, BackpressureConfig cfg)
      : sink_(std::move(sink)), cfg_(cfg), queue_(cfg_.max_queue_msgs, cfg_.max_queue_bytes) {
    validate_config(cfg_);
  }

  ~SinkIngress() {
    // Handle orderly shutdown to avoid std::terminate on joinable threads
    queue_.close();
    for (auto &w : writers_) {
      if (w.joinable()) w.join();
    }
  }

  void start() {
    if (start_called_) throw std::logic_error("SinkIngress::start called twice");
    start_called_      = true;
    const auto threads = cfg_.writer_threads;
    active_writers_.store(threads, std::memory_order_relaxed);
    writers_.reserve(threads);
    for (std::size_t i = 0; i < threads; ++i) {
      writers_.emplace_back([this] { this->run(); });
    }
  }

  void run() {
    std::vector<QueueNode> batch;
    batch.reserve(cfg_.max_batch_msgs);
    try {
      for (;;) {
        batch.clear();
        if (!queue_.pop_batch(batch, cfg_.max_batch_msgs, cfg_.max_batch_bytes)) break;
        for (const auto &n : batch) {
          EmissionView v{
              n.channel_id,
              std::string_view{n.payload->data(), n.payload->size()}
          };
          sink_->write(v);
          ++written_msgs_;
          written_bytes_ += n.size;
        }
      }
    } catch (const std::exception &ex) {
      record_failure(ex.what());
    } catch (...) {
      record_failure("unknown error in sink writer");
    }
    finalize_writer();
  }

  void close_queue() { queue_.close(); }

  bool join_until(std::chrono::steady_clock::time_point deadline) {
    // Always acquire mutex before checking finished_ to avoid missing notifications
    std::unique_lock<std::mutex> lk(finished_m_);
    if (!finished_cv_.wait_until(lk, deadline, [&] { return finished_.load(std::memory_order_acquire); })) {
      return false; // deadline reached
    }
    lk.unlock();

    for (auto &w : writers_) {
      if (w.joinable()) w.join();
    }
    return true;
  }

  bool offer(uint32_t ch_id, std::shared_ptr<const std::string> buf, std::size_t sz) {
    QueueNode n{ch_id, std::move(buf), sz};
    const bool ok = queue_.offer(std::move(n), cfg_, &stall_ns_accum_, &drops_counter_);
    if (ok) {
      enq_msgs_.fetch_add(1, std::memory_order_relaxed);
      enq_bytes_.fetch_add(sz, std::memory_order_relaxed);
    }
    return ok;
  }

  // ingress state
  std::shared_ptr<IDynamicSink> sink_;
  // Config is captured at construction time and remains fixed for the sink's lifetime.
  // Subsequent default changes do not affect existing ingresses.
  BackpressureConfig cfg_{};
  BoundedQueue queue_;
  std::vector<std::thread> writers_;

  // clang-format off
  // metrics/state
  alignas(64) std::atomic<uint64_t> enq_msgs_      {0};
  alignas(64) std::atomic<uint64_t> enq_bytes_     {0};
  alignas(64) std::atomic<uint64_t> written_msgs_  {0};
  alignas(64) std::atomic<uint64_t> written_bytes_ {0};
  alignas(64) std::atomic<uint64_t> stall_ns_accum_{0};
  alignas(64) std::atomic<uint64_t> drops_counter_ {0};
  alignas(64) std::atomic<bool>     finished_  {false};
  alignas(64) std::atomic<bool>     failed_    {false};
  alignas(64) std::atomic<std::size_t> active_writers_{0};
  // clang-format on
  // error_ is written by the writer thread and only read after join_until() observes finished_ true
  // and joins the thread.
  std::string error_;

  std::condition_variable finished_cv_;
  mutable std::mutex finished_m_;

  bool start_called_ = false;

  void finalize_writer() {
    if (active_writers_.fetch_sub(1, std::memory_order_acq_rel) == 1) {
      try {
        sink_->flush();
      } catch (const std::exception &ex) {
        record_failure(ex.what());
      } catch (...) {
        record_failure("unknown error in sink writer");
      }
      try {
        sink_->close();
      } catch (const std::exception &ex) {
        record_failure(ex.what());
      } catch (...) {
        record_failure("unknown error in sink writer");
      }
      {
        std::lock_guard<std::mutex> lk(finished_m_);
        finished_.store(true, std::memory_order_release);
      }
      finished_cv_.notify_all();
    }
  }

  void record_failure(std::string msg) {
    queue_.close();
    bool expected = false;
    if (failed_.compare_exchange_strong(expected, true, std::memory_order_relaxed)) {
      error_ = std::move(msg);
    }
  }

  // For observability
  uint64_t drops() const { return drops_counter_.load(std::memory_order_relaxed); }
  uint64_t stalled_ns() const { return stall_ns_accum_.load(std::memory_order_relaxed); }
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_BACKPRESSURE_HPP
