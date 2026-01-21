#ifndef LAHUTA_PIPELINE_DYNAMIC_PROGRESS_OBSERVER_HPP
#define LAHUTA_PIPELINE_DYNAMIC_PROGRESS_OBSERVER_HPP

#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <mutex>

#include <spdlog/sinks/sink.h>

#include "pipeline/dynamic/run_observer.hpp"

// clang-format off
namespace lahuta::pipeline::dynamic {

struct ProgressObserverConfig {
  std::chrono::milliseconds interval{50};
};

class ProgressRunObserver final : public IRunObserver {
  friend class ProgressAwareSink;
public:
  explicit ProgressRunObserver(ProgressObserverConfig config = {})
      : output_(std::cerr), enabled_(config.interval.count() > 0) {
    auto interval = config.interval;
    if (interval < std::chrono::milliseconds{0}) interval = std::chrono::milliseconds{0};
    interval_ns_ = std::chrono::duration_cast<std::chrono::nanoseconds>(interval).count();
    const auto now = Clock::now();
    next_report_ns_.store(to_ns(now) + interval_ns_, std::memory_order_relaxed);
  }

  void on_item_end(std::size_t, const PipelineItem &) override {
    if (!enabled_) return;
    const auto count = processed_.fetch_add(1, std::memory_order_relaxed) + 1;
    maybe_report(count);
  }

  void on_item_begin(std::size_t, const PipelineItem &) override {
    if (!enabled_) return;
    bool expected = false;
    if (!started_.compare_exchange_strong(expected, true, std::memory_order_relaxed)) return;
    const auto now = Clock::now();
    next_report_ns_.store(to_ns(now) + interval_ns_, std::memory_order_relaxed);
    report(0);
  }

  void finish() {
    if (!enabled_) return;
    if (finished_.exchange(true, std::memory_order_relaxed)) return;
    clear_line();
  }

  std::uint64_t processed() const noexcept {
    return processed_.load(std::memory_order_relaxed);
  }

private:
  using Clock = std::chrono::steady_clock;
  static constexpr const char Label[] = "processed";
  static constexpr std::size_t LabelLength = sizeof(Label) - 1;

  static std::int64_t to_ns(Clock::time_point tp) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(tp.time_since_epoch()).count();
  }

  static std::size_t digit_count(std::uint64_t value) {
    std::size_t digits = 1;
    while (value >= 10) {
      value /= 10;
      ++digits;
    }
    return digits;
  }

  void maybe_report(std::uint64_t count) {
    if (interval_ns_ <= 0) return;
    const auto now = Clock::now();
    const auto now_ns = to_ns(now);
    auto next_ns = next_report_ns_.load(std::memory_order_relaxed);
    if (now_ns < next_ns) return;
    if (!next_report_ns_.compare_exchange_strong(next_ns, now_ns + interval_ns_, std::memory_order_relaxed)) return;
    report(count);
  }

  void report(std::uint64_t count) {
    const auto line_len = LabelLength + 1 + digit_count(count);
    std::lock_guard<std::mutex> lock(output_mutex_);
    output_ << '\r' << Label << '=' << count << std::flush;
    last_count_ = count;
    last_len_ = line_len;
  }

  void clear_line() {
    auto lock = lock_output();
    clear_line_locked(true);
  }

  std::unique_lock<std::mutex> lock_output() {
    return std::unique_lock<std::mutex>(output_mutex_);
  }

  void clear_line_locked(bool reset) {
    if (last_len_ == 0) return;
    output_ << '\r';
    for (std::size_t i = 0; i < last_len_; ++i) output_ << ' ';
    output_ << '\r' << std::flush;
    if (reset) last_len_ = 0;
  }

  void reprint_locked() {
    if (last_len_ == 0) return;
    output_ << '\r' << Label << '=' << last_count_ << std::flush;
  }

  std::atomic<std::uint64_t> processed_{0};
  std::atomic<std::int64_t> next_report_ns_{0};
  std::int64_t interval_ns_{0};
  std::ostream& output_;
  const bool enabled_{true};
  std::atomic<bool> finished_{false};
  std::atomic<bool> started_{false};
  std::mutex output_mutex_;
  std::size_t last_len_{0};
  std::uint64_t last_count_{0};
};

class ProgressAwareSink : public spdlog::sinks::sink {
public:
  ProgressAwareSink(std::shared_ptr<spdlog::sinks::sink> sink, std::shared_ptr<ProgressRunObserver> observer)
      : sink_(std::move(sink)), observer_(std::move(observer)) {}

  void log(const spdlog::details::log_msg &msg) override {
    auto observer = observer_.lock();
    if (!observer) {
      sink_->log(msg);
      return;
    }
    auto lock = observer->lock_output();
    observer->clear_line_locked(false);
    sink_->log(msg);
    observer->reprint_locked();
  }

  void flush() override { sink_->flush(); }
  void set_pattern(const std::string &pattern) override { sink_->set_pattern(pattern); }
  void set_formatter(std::unique_ptr<spdlog::formatter> sink_formatter) override {
    sink_->set_formatter(std::move(sink_formatter));
  }

  std::shared_ptr<spdlog::sinks::sink> wrapped_sink() const { return sink_; }

private:
  std::shared_ptr<spdlog::sinks::sink> sink_;
  std::weak_ptr<ProgressRunObserver> observer_;
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_PROGRESS_OBSERVER_HPP
