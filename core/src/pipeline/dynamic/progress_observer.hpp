#ifndef LAHUTA_PIPELINE_DYNAMIC_PROGRESS_OBSERVER_HPP
#define LAHUTA_PIPELINE_DYNAMIC_PROGRESS_OBSERVER_HPP

#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <mutex>
#include <optional>
#include <string>

#include <spdlog/fmt/fmt.h>
#include <spdlog/sinks/sink.h>

#include "pipeline/dynamic/run_observer.hpp"

// clang-format off
namespace lahuta::pipeline::dynamic {

struct ProgressObserverConfig {
  std::chrono::milliseconds interval{50};
  std::string label{"progress"};
  std::optional<std::size_t> total_items;
};

class ProgressRunObserver final : public IRunObserver {
  friend class ProgressAwareSink;
public:
  explicit ProgressRunObserver(ProgressObserverConfig config = {})
      : output_(std::cerr),
        label_(config.label.empty() ? "progress" : std::move(config.label)),
        total_items_(config.total_items),
        enabled_(config.interval.count() > 0) {
    auto interval = config.interval;
    if (interval < std::chrono::milliseconds{0}) interval = std::chrono::milliseconds{0};
    interval_ns_ = std::chrono::duration_cast<std::chrono::nanoseconds>(interval).count();
    const auto now = Clock::now();
    next_report_ns_.store(to_ns(now) + interval_ns_, std::memory_order_relaxed);
  }

  void on_item_end(std::size_t, const PipelineItem &) override {
    if (!enabled_) return;
    inflight_.fetch_sub(1, std::memory_order_relaxed);
    const auto count = processed_.fetch_add(1, std::memory_order_relaxed) + 1;
    maybe_report(count);
  }

  void on_item_begin(std::size_t, const PipelineItem &) override {
    if (!enabled_) return;
    inflight_.fetch_add(1, std::memory_order_relaxed);
    bool expected = false;
    if (!started_.compare_exchange_strong(expected, true, std::memory_order_relaxed)) return;
    const auto now = Clock::now();
    next_report_ns_.store(to_ns(now) + interval_ns_, std::memory_order_relaxed);
    {
      std::lock_guard<std::mutex> lock(output_mutex_);
      start_time_ = now;
      last_rate_sample_time_ = now;
      last_rate_sample_count_ = 0;
      last_rate_items_per_sec_ = 0.0;
      report_locked(0, now);
    }
  }

  void on_item_skipped(std::size_t, const PipelineItem &, std::string_view) override {
    if (!enabled_) return;
    skipped_.fetch_add(1, std::memory_order_relaxed);
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
  static constexpr double RateSampleSeconds = 0.5;

  static std::int64_t to_ns(Clock::time_point tp) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(tp.time_since_epoch()).count();
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
    std::lock_guard<std::mutex> lock(output_mutex_);
    report_locked(count, Clock::now());
  }

  void report_locked(std::uint64_t count, Clock::time_point now) {
    const auto skipped = skipped_.load(std::memory_order_relaxed);
    const auto ok = count >= skipped ? count - skipped : 0;
    const auto inflight = inflight_.load(std::memory_order_relaxed);
    const double elapsed = std::chrono::duration<double>(now - start_time_).count();
    const double avg_rate = elapsed > 0.0 ? static_cast<double>(count) / elapsed : 0.0;

    const auto rate_dt = std::chrono::duration<double>(now - last_rate_sample_time_).count();
    if (rate_dt >= RateSampleSeconds) {
      if (count >= last_rate_sample_count_) {
        last_rate_items_per_sec_ = static_cast<double>(count - last_rate_sample_count_) / rate_dt;
      }
      last_rate_sample_count_ = count;
      last_rate_sample_time_ = now;
    }
    const double rate = last_rate_items_per_sec_ > 0.0 ? last_rate_items_per_sec_ : avg_rate;

    buffer_.clear();
    fmt::format_to(std::back_inserter(buffer_), "{}: done={} ok={} skip={} inflight={} rate={:.2f}/s elapsed=",
                   label_, count, ok, skipped, inflight, rate);
    append_duration(buffer_, elapsed);

    if (total_items_ && *total_items_ > 0) {
      const auto total = *total_items_;
      const double pct = static_cast<double>(count) * 100.0 / static_cast<double>(total);
      fmt::format_to(std::back_inserter(buffer_), " total={} ({:.1f}%)", total, pct);
      if (rate > 0.0 && count < total) {
        const double eta = static_cast<double>(total - count) / rate;
        fmt::format_to(std::back_inserter(buffer_), " eta=");
        append_duration(buffer_, eta);
      }
    }

    last_line_.assign(buffer_.data(), buffer_.size());
    output_ << '\r' << last_line_ << std::flush;
    last_len_ = last_line_.size();
  }

  static void append_duration(fmt::memory_buffer &buf, double seconds) {
    if (seconds < 0.0) seconds = 0.0;
    const auto total_ms = static_cast<std::uint64_t>(seconds * 1000.0);
    const auto ms = total_ms % 1000;
    const auto total_s = total_ms / 1000;
    const auto s = total_s % 60;
    const auto m = (total_s / 60) % 60;
    const auto h = total_s / 3600;
    if (h > 0) {
      fmt::format_to(std::back_inserter(buf), "{:02}:{:02}:{:02}", h, m, s);
    } else {
      fmt::format_to(std::back_inserter(buf), "{:02}:{:02}.{:01}", m, s, ms / 100);
    }
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
    if (reset) {
      last_len_ = 0;
      last_line_.clear();
    }
  }

  void reprint_locked() {
    if (last_len_ == 0 || last_line_.empty()) return;
    output_ << '\r' << last_line_ << std::flush;
  }

  std::atomic<std::uint64_t> processed_{0};
  std::atomic<std::uint64_t> skipped_{0};
  std::atomic<std::int64_t> inflight_{0};
  std::atomic<std::int64_t> next_report_ns_{0};
  std::int64_t interval_ns_{0};
  std::ostream& output_;
  std::string label_;
  std::optional<std::size_t> total_items_;
  const bool enabled_{true};
  std::atomic<bool> finished_{false};
  std::atomic<bool> started_{false};
  std::mutex output_mutex_;
  Clock::time_point start_time_{};
  Clock::time_point last_rate_sample_time_{};
  std::uint64_t last_rate_sample_count_{0};
  double last_rate_items_per_sec_{0.0};
  fmt::memory_buffer buffer_;
  std::size_t last_len_{0};
  std::string last_line_;
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
