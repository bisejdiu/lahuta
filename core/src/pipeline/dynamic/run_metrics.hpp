#ifndef LAHUTA_PIPELINE_DYNAMIC_RUN_METRICS_HPP
#define LAHUTA_PIPELINE_DYNAMIC_RUN_METRICS_HPP

#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <utility>
#include <vector>

namespace lahuta::pipeline::dynamic {

struct StageMetricsSnapshot {
  std::int64_t ingest_ns  = 0;
  std::int64_t prepare_ns = 0;
  std::int64_t setup_ns   = 0;
  std::int64_t compute_ns = 0;
  std::int64_t flush_ns   = 0;
  std::size_t  items_total   = 0;
  std::size_t  items_skipped = 0;
};

class StageRunMetrics {
  struct Slot;
public:
  using Clock = std::chrono::steady_clock;

  class ThreadHandle {
  public:
    ThreadHandle() = default;
    ThreadHandle(const ThreadHandle&) = delete;
    ThreadHandle& operator=(const ThreadHandle&) = delete;

    ThreadHandle(ThreadHandle&& other) noexcept { *this = std::move(other); }
    ThreadHandle& operator=(ThreadHandle&& other) noexcept {
      if (this != &other) {
        reset();
        owner_ = other.owner_;
        slot_  = other.slot_;
        other.owner_ = nullptr;
        other.slot_  = nullptr;
      }
      return *this;
    }

    ~ThreadHandle() { reset(); }

    bool attached() const noexcept {
      return owner_ != nullptr && slot_ != nullptr;
    }

    void reset() noexcept {
      owner_ = nullptr;
      slot_  = nullptr;
    }

  private:
    friend class StageRunMetrics;

    StageRunMetrics* owner_{nullptr};
    Slot*            slot_{nullptr};
  };

  StageRunMetrics() = default;
  StageRunMetrics(const StageRunMetrics&) = delete;
  StageRunMetrics& operator=(const StageRunMetrics&) = delete;

  void ensure(ThreadHandle& handle) {
    if (handle.owner_ != this || handle.slot_ == nullptr) {
      allocate_slot(handle);
    }
  }

  void add_ingest(ThreadHandle& handle, std::chrono::nanoseconds d) {
    auto* slot = slot_for(handle);
    slot->ingest_ns += d.count();
  }

  void add_prepare(ThreadHandle& handle, std::chrono::nanoseconds d) {
    auto* slot = slot_for(handle);
    slot->prepare_ns += d.count();
  }

  void add_setup(ThreadHandle& handle, std::chrono::nanoseconds d) {
    auto* slot = slot_for(handle);
    slot->setup_ns += d.count();
  }

  void add_compute(ThreadHandle& handle, std::chrono::nanoseconds d) {
    auto* slot = slot_for(handle);
    slot->compute_ns += d.count();
  }

  void add_flush(std::chrono::nanoseconds d) {
    flush_ns_.fetch_add(d.count(), std::memory_order_relaxed);
  }

  void inc_items_total(ThreadHandle& handle) {
    auto* slot = slot_for(handle);
    ++slot->items_total;
  }

  void inc_items_skipped(ThreadHandle& handle) {
    auto* slot = slot_for(handle);
    ++slot->items_skipped;
  }

  StageMetricsSnapshot snapshot() const {
    StageMetricsSnapshot out;
    std::scoped_lock lk(mutex_);
    for (const auto& ptr : slots_) {
      if (!ptr) continue;
      out.ingest_ns     += ptr->ingest_ns;
      out.prepare_ns    += ptr->prepare_ns;
      out.setup_ns      += ptr->setup_ns;
      out.compute_ns    += ptr->compute_ns;
      out.items_total   += ptr->items_total;
      out.items_skipped += ptr->items_skipped;
    }
    out.flush_ns = flush_ns_.load(std::memory_order_relaxed);
    return out;
  }

private:
  struct Slot {
    std::int64_t ingest_ns  = 0;
    std::int64_t prepare_ns = 0;
    std::int64_t setup_ns   = 0;
    std::int64_t compute_ns = 0;
    std::size_t  items_total   = 0;
    std::size_t  items_skipped = 0;
  };

  Slot* slot_for(ThreadHandle& handle) {
    ensure(handle);
    return handle.slot_;
  }

  void allocate_slot(ThreadHandle& handle) {
    auto slot_ptr = std::make_unique<Slot>();
    Slot* raw = slot_ptr.get();
    {
      std::scoped_lock lk(mutex_);
      slots_.push_back(std::move(slot_ptr));
    }
    handle.owner_ = this;
    handle.slot_  = raw;
  }

  mutable std::mutex mutex_;
  std::vector<std::unique_ptr<Slot>> slots_;
  std::atomic<std::int64_t> flush_ns_{0};
};

class NullStageRunMetrics {
public:
  using Clock = std::chrono::steady_clock;

  class ThreadHandle {
  public:
    void reset() noexcept {}
  };

  void ensure(ThreadHandle&) noexcept {}
  void add_ingest(ThreadHandle&, std::chrono::nanoseconds) noexcept {}
  void add_prepare(ThreadHandle&, std::chrono::nanoseconds) noexcept {}
  void add_setup(ThreadHandle&, std::chrono::nanoseconds) noexcept {}
  void add_compute(ThreadHandle&, std::chrono::nanoseconds) noexcept {}
  void add_flush(std::chrono::nanoseconds) noexcept {}
  void inc_items_total(ThreadHandle&) noexcept {}
  void inc_items_skipped(ThreadHandle&) noexcept {}

  StageMetricsSnapshot snapshot() const noexcept { return {}; }
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_RUN_METRICS_HPP
