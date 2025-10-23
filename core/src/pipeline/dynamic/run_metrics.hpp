#ifndef LAHUTA_PIPELINE_DYNAMIC_RUN_METRICS_HPP
#define LAHUTA_PIPELINE_DYNAMIC_RUN_METRICS_HPP

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <utility>
#include <vector>

// clang-format off
namespace lahuta::pipeline::dynamic {

struct StageMetricsSnapshot {
  std::int64_t ingest_ns  = 0;
  std::int64_t prepare_ns = 0;
  std::int64_t setup_ns   = 0;
  std::int64_t compute_ns = 0;
  std::int64_t flush_ns   = 0;
  std::size_t  items_total   = 0;
  std::size_t  items_skipped = 0;
  std::uint64_t permit_wait_ns_total = 0;
  std::uint64_t permit_wait_samples  = 0;
  std::int64_t  permit_wait_ns_min   = 0;
  std::int64_t  permit_wait_ns_max   = 0;
  std::size_t   inflight_peak        = 0;
  std::uint64_t inflight_sum         = 0;
  std::uint64_t inflight_samples     = 0;
  std::vector<std::int64_t> stage_setup_ns;
  std::vector<std::int64_t> stage_compute_ns;
};

class StageRunMetrics {
  struct Slot;

public:
  using Clock = std::chrono::steady_clock;

  class ThreadHandle {
  public:
    ThreadHandle() = default;
    ThreadHandle(const ThreadHandle &) = delete;
    ThreadHandle &operator=(const ThreadHandle &) = delete;

    ThreadHandle(ThreadHandle &&other) noexcept { *this = std::move(other); }
    ThreadHandle &operator=(ThreadHandle &&other) noexcept {
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

    bool attached() const noexcept { return owner_ != nullptr && slot_ != nullptr; }

    void reset() noexcept {
      owner_ = nullptr;
      slot_  = nullptr;
    }

  private:
    friend class StageRunMetrics;

    StageRunMetrics *owner_{nullptr};
    Slot *slot_{nullptr};
  };

  StageRunMetrics() = default;
  explicit StageRunMetrics(bool enable_stage_breakdown) : stage_breakdown_enabled_(enable_stage_breakdown) {}
  StageRunMetrics(const StageRunMetrics &) = delete;
  StageRunMetrics &operator=(const StageRunMetrics &) = delete;

  void set_stage_breakdown_enabled(bool enabled) {
    stage_breakdown_enabled_.store(enabled, std::memory_order_relaxed);
  }

  bool stage_breakdown_enabled() const noexcept {
    return stage_breakdown_enabled_.load(std::memory_order_relaxed);
  }

  void configure_stage_breakdown(std::size_t count) {
    if (!stage_breakdown_enabled()) return;
    std::scoped_lock lk(mutex_);
    stage_count_ = count;
    for (auto &ptr : slots_) {
      if (ptr) ptr->ensure_stage_capacity(stage_count_);
    }
  }

  void ensure(ThreadHandle &handle) {
    if (handle.owner_ != this || handle.slot_ == nullptr) {
      allocate_slot(handle);
    }
  }

  void add_ingest(ThreadHandle &handle, std::chrono::nanoseconds d) {
    auto *slot = slot_for(handle);
    slot->ingest_ns += d.count();
  }

  void add_prepare(ThreadHandle &handle, std::chrono::nanoseconds d) {
    auto *slot = slot_for(handle);
    slot->prepare_ns += d.count();
  }

  void add_setup(ThreadHandle &handle, std::chrono::nanoseconds d) {
    auto *slot = slot_for(handle);
    slot->setup_ns += d.count();
  }

  void add_compute(ThreadHandle &handle, std::chrono::nanoseconds d) {
    auto *slot = slot_for(handle);
    slot->compute_ns += d.count();
  }

  void add_flush(std::chrono::nanoseconds d) { flush_ns_.fetch_add(d.count(), std::memory_order_relaxed); }

  void inc_items_total(ThreadHandle &handle) {
    auto *slot = slot_for(handle);
    ++slot->items_total;
  }

  void inc_items_skipped(ThreadHandle &handle) {
    auto *slot = slot_for(handle);
    ++slot->items_skipped;
  }

  void add_permit_wait(ThreadHandle &handle, std::chrono::nanoseconds d) {
    if (d.count() < 0) return;
    auto *slot = slot_for(handle);
    const auto ns = d.count();
    slot->permit_wait_ns_total += static_cast<std::uint64_t>(ns);
    ++slot->permit_wait_samples;
    if (ns < slot->permit_wait_ns_min) slot->permit_wait_ns_min = ns;
    if (ns > slot->permit_wait_ns_max) slot->permit_wait_ns_max = ns;
  }

  void on_item_inflight_enter() {
    const auto current = inflight_current_.fetch_add(1, std::memory_order_relaxed) + 1;
    update_peak(current);
    inflight_sum_.fetch_add(current, std::memory_order_relaxed);
    inflight_samples_.fetch_add(1, std::memory_order_relaxed);
  }

  void on_item_inflight_exit() {
    const auto previous = inflight_current_.fetch_sub(1, std::memory_order_relaxed);
    const auto current = previous > 0 ? previous - 1 : 0;
    inflight_sum_.fetch_add(current, std::memory_order_relaxed);
    inflight_samples_.fetch_add(1, std::memory_order_relaxed);
  }

  void add_stage_setup(ThreadHandle &handle, std::size_t stage_index, std::chrono::nanoseconds d) {
    if (!stage_breakdown_enabled()) return;
    auto *slot = slot_for(handle);
    if (stage_index >= slot->stage_setup_ns.size()) return;
    slot->stage_setup_ns[stage_index] += d.count();
  }

  void add_stage_compute(ThreadHandle &handle, std::size_t stage_index, std::chrono::nanoseconds d) {
    if (!stage_breakdown_enabled()) return;
    auto *slot = slot_for(handle);
    if (stage_index >= slot->stage_compute_ns.size()) return;
    slot->stage_compute_ns[stage_index] += d.count();
  }

  StageMetricsSnapshot snapshot() const {
    StageMetricsSnapshot out;
    std::scoped_lock lk(mutex_);
    if (stage_breakdown_enabled()) {
      out.stage_setup_ns.assign(stage_count_, 0);
      out.stage_compute_ns.assign(stage_count_, 0);
    }
    std::int64_t global_min_wait = std::numeric_limits<std::int64_t>::max();
    std::int64_t global_max_wait = 0;
    std::uint64_t total_wait_ns  = 0;
    std::uint64_t wait_samples   = 0;
    for (const auto &ptr : slots_) {
      if (!ptr) continue;
      out.ingest_ns   += ptr->ingest_ns;
      out.prepare_ns  += ptr->prepare_ns;
      out.setup_ns    += ptr->setup_ns;
      out.compute_ns  += ptr->compute_ns;
      out.items_total += ptr->items_total;
      out.items_skipped += ptr->items_skipped;
      if (ptr->permit_wait_samples > 0) {
        total_wait_ns += ptr->permit_wait_ns_total;
        wait_samples  += ptr->permit_wait_samples;
        if (ptr->permit_wait_ns_min < global_min_wait) global_min_wait = ptr->permit_wait_ns_min;
        if (ptr->permit_wait_ns_max > global_max_wait) global_max_wait = ptr->permit_wait_ns_max;
      }
      if (stage_breakdown_enabled()) {
        const auto n = std::min(stage_count_, ptr->stage_setup_ns.size());
        for (std::size_t i = 0; i < n; ++i) {
          out.stage_setup_ns[i]   += ptr->stage_setup_ns[i];
          out.stage_compute_ns[i] += ptr->stage_compute_ns[i];
        }
      }
    }
    out.permit_wait_ns_total = total_wait_ns;
    out.permit_wait_samples  = wait_samples;
    if (wait_samples > 0) {
      out.permit_wait_ns_min = global_min_wait;
      out.permit_wait_ns_max = global_max_wait;
    } else {
      out.permit_wait_ns_min = 0;
      out.permit_wait_ns_max = 0;
    }
    out.flush_ns = flush_ns_.load(std::memory_order_relaxed);
    out.inflight_peak = inflight_peak_.load(std::memory_order_relaxed);
    out.inflight_sum  = inflight_sum_.load(std::memory_order_relaxed);
    out.inflight_samples = inflight_samples_.load(std::memory_order_relaxed);
    return out;
  }

private:
  struct Slot {
    std::int64_t ingest_ns  = 0;
    std::int64_t prepare_ns = 0;
    std::int64_t setup_ns   = 0;
    std::int64_t compute_ns = 0;
    std::size_t items_total = 0;
    std::size_t items_skipped = 0;
    std::uint64_t permit_wait_ns_total = 0;
    std::uint64_t permit_wait_samples  = 0;
    std::int64_t permit_wait_ns_min = std::numeric_limits<std::int64_t>::max();
    std::int64_t permit_wait_ns_max = 0;
    std::vector<std::int64_t> stage_setup_ns;
    std::vector<std::int64_t> stage_compute_ns;

    void ensure_stage_capacity(std::size_t count) {
      stage_setup_ns.assign(count, 0);
      stage_compute_ns.assign(count, 0);
    }
  };

  Slot *slot_for(ThreadHandle &handle) {
    ensure(handle);
    return handle.slot_;
  }

  void allocate_slot(ThreadHandle &handle) {
    auto slot_ptr = std::make_unique<Slot>();
    Slot *raw = slot_ptr.get();
    {
      std::scoped_lock lk(mutex_);
      if (stage_breakdown_enabled()) {
        slot_ptr->ensure_stage_capacity(stage_count_);
      }
      slots_.push_back(std::move(slot_ptr));
    }
    handle.owner_ = this;
    handle.slot_ = raw;
  }

  void update_peak(std::size_t value) {
    std::size_t expected = inflight_peak_.load(std::memory_order_relaxed);
    while (value > expected && !inflight_peak_.compare_exchange_weak(expected, value, std::memory_order_relaxed)) {
      // expected is updated with current value on failure
    }
  }

  mutable std::mutex mutex_;
  std::vector<std::unique_ptr<Slot>> slots_;
  std::atomic<std::int64_t> flush_ns_{0};
  std::atomic<std::size_t>  inflight_current_{0};
  std::atomic<std::size_t>  inflight_peak_{0};
  std::atomic<std::uint64_t> inflight_sum_{0};
  std::atomic<std::uint64_t> inflight_samples_{0};
  std::atomic<bool> stage_breakdown_enabled_{false};
  std::size_t stage_count_{0};
};

class NullStageRunMetrics {
public:
  using Clock = std::chrono::steady_clock;

  class ThreadHandle {
  public:
    void reset() noexcept {}
  };

  void set_stage_breakdown_enabled(bool) noexcept {}
  bool stage_breakdown_enabled() const noexcept { return false; }
  void configure_stage_breakdown(std::size_t) noexcept {}

  void ensure(ThreadHandle &) noexcept {}
  void add_ingest(ThreadHandle &, std::chrono::nanoseconds) noexcept {}
  void add_prepare(ThreadHandle &, std::chrono::nanoseconds) noexcept {}
  void add_setup(ThreadHandle &, std::chrono::nanoseconds) noexcept {}
  void add_compute(ThreadHandle &, std::chrono::nanoseconds) noexcept {}
  void add_flush(std::chrono::nanoseconds) noexcept {}
  void inc_items_total(ThreadHandle &) noexcept {}
  void inc_items_skipped(ThreadHandle &) noexcept {}
  void add_permit_wait(ThreadHandle &, std::chrono::nanoseconds) noexcept {}
  void on_item_inflight_enter() noexcept {}
  void on_item_inflight_exit() noexcept {}
  void add_stage_setup(ThreadHandle &, std::size_t, std::chrono::nanoseconds) noexcept {}
  void add_stage_compute(ThreadHandle &, std::size_t, std::chrono::nanoseconds) noexcept {}

  StageMetricsSnapshot snapshot() const noexcept { return {}; }
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_RUN_METRICS_HPP
