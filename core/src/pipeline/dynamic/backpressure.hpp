#ifndef LAHUTA_PIPELINE_DYNAMIC_BACKPRESSURE_HPP
#define LAHUTA_PIPELINE_DYNAMIC_BACKPRESSURE_HPP

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

#include "pipeline/dynamic/sink_iface.hpp"
#include "pipeline/dynamic/types.hpp"

// clang-format off
namespace lahuta::pipeline::dynamic {

// Backpressure policy when a per-sink queue is full
enum class OnFull { Block, DropLatest, DropOldest };

struct BackpressureConfig {
  std::size_t max_queue_msgs  = 65536;
  std::size_t max_queue_bytes = 256 * 1024 * 1024;
  std::size_t max_batch_msgs  = 256;
  std::size_t max_batch_bytes = 1 * 1024 * 1024;
  std::chrono::milliseconds offer_timeout{500};
  OnFull on_full = OnFull::Block;
  bool required = true; // failure aborts pipeline
};

// Node stored in per-sink queue
struct QueueNode {
  uint32_t channel_id = 0;
  std::shared_ptr<const std::string> payload;
  std::size_t size = 0;
};

class BoundedQueue {
public:
  BoundedQueue(std::size_t max_msgs, std::size_t max_bytes) : max_msgs_(max_msgs), max_bytes_(max_bytes) {}

  bool offer(QueueNode n,
             const BackpressureConfig& cfg,
             std::atomic<uint64_t>* stall_ns_accum,
             std::atomic<uint64_t>* drops_counter) {
    const auto need = n.size;
    std::unique_lock<std::mutex> lk(m_);
    for (;;) {
      if (closed_) return false;
      if ((q_.size() < max_msgs_) && (bytes_ + need <= max_bytes_)) {
        q_.push_back(std::move(n));
        bytes_ += need;
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
          cv_not_full_.wait_for(lk, cfg.offer_timeout, [&]{
            return closed_ || ((q_.size() < max_msgs_) && (bytes_ + need <= max_bytes_));
          });
          if (closed_) return false;
          auto t1 = std::chrono::steady_clock::now();
          if (stall_ns_accum) *stall_ns_accum += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
          continue;
        }
      }
    }
  }

  bool pop_batch(std::vector<QueueNode>& out, std::size_t max_msgs, std::size_t max_bytes) {
    std::unique_lock<std::mutex> lk(m_);
    cv_not_empty_.wait(lk, [&]{ return closed_ || !q_.empty(); });
    if (q_.empty() && closed_) return false;
    std::size_t taken_msgs = 0, taken_bytes = 0;
    while (!q_.empty() && taken_msgs < max_msgs) {
      const auto& front = q_.front();
      if (taken_bytes + front.size > max_bytes) break;
      taken_bytes += front.size;
      ++taken_msgs;
      out.push_back(front);
      bytes_ -= front.size;
      q_.pop_front();
    }
    lk.unlock();
    cv_not_full_.notify_all();
    return !out.empty();
  }

  void close() {
    std::lock_guard<std::mutex> lk(m_);
    closed_ = true;
    cv_not_empty_.notify_all();
    cv_not_full_ .notify_all();
  }

  std::size_t size_msgs() const {
    std::lock_guard<std::mutex> lk(m_);
    return q_.size();
  }

  std::size_t size_bytes() const {
    std::lock_guard<std::mutex> lk(m_);
    return bytes_;
  }

private:
  std::size_t max_msgs_ = 0;
  std::size_t max_bytes_ = 0;
  std::deque<QueueNode> q_;
  std::size_t bytes_ = 0;
  mutable std::mutex m_;
  std::condition_variable cv_not_full_;
  std::condition_variable cv_not_empty_;
  bool closed_ = false;
};

struct SinkIngress {
  explicit SinkIngress(std::shared_ptr<IDynamicSink> sink, BackpressureConfig cfg)
      : sink_(std::move(sink)), cfg_(cfg), queue_(cfg_.max_queue_msgs, cfg_.max_queue_bytes) {}

  ~SinkIngress() {
    // Handle orderly shutdown to avoid std::terminate on joinable threads
    queue_.close();
    if (writer_.joinable()) writer_.join();
  }

  void start() {
    writer_ = std::thread([this]{ this->run(); });
  }

  void run() {
    std::vector<QueueNode> batch;
    try {
      for (;;) {
        batch.clear();
        if (!queue_.pop_batch(batch, cfg_.max_batch_msgs, cfg_.max_batch_bytes)) break;
        for (const auto& n : batch) {
          EmissionView v{n.channel_id, std::string_view{n.payload->data(), n.payload->size()}};
          sink_->write(v);
          ++written_msgs_;
          written_bytes_ += n.size;
        }
      }
      sink_->flush();
      sink_->close();
    } catch (const std::exception& ex) {
      failed_ = true;
      error_ = ex.what();
    } catch (...) {
      failed_ = true;
      error_ = "unknown error in sink writer";
    }
    finished_ = true;
  }

  void close_queue() { queue_.close(); }

  bool join_until(std::chrono::steady_clock::time_point deadline) {
    // Poll finished_ to avoid blocking join beyond deadline
    while (!finished_) {
      if (std::chrono::steady_clock::now() >= deadline) return false;
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    if (writer_.joinable()) writer_.join();
    return true;
  }

  bool offer(uint32_t ch_id, std::shared_ptr<const std::string> buf, std::size_t sz) {
    QueueNode n{ch_id, std::move(buf), sz};
    return queue_.offer(std::move(n), cfg_, &stall_ns_accum_, &drops_counter_);
  }

  // ingress state
  std::shared_ptr<IDynamicSink> sink_;
  BackpressureConfig cfg_{};
  BoundedQueue queue_;
  std::thread writer_;

  // metrics/state
  std::atomic<uint64_t> enq_msgs_      {0};
  std::atomic<uint64_t> enq_bytes_     {0};
  std::atomic<uint64_t> written_msgs_  {0};
  std::atomic<uint64_t> written_bytes_ {0};
  std::atomic<uint64_t> stall_ns_accum_{0};
  std::atomic<uint64_t> drops_counter_ {0};
  std::atomic<bool>     finished_  {false};
  std::atomic<bool>     failed_    {false};
  std::string           error_;
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_BACKPRESSURE_HPP
