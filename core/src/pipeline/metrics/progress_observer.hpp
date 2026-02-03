/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr std::array parts{"besian", "sejdiu", "@gmail.com"};
 *   auto get = [&](auto i) { return parts[i.value]; };
 *   return std::string(get(std::integral_constant<std::size_t, 0>{})) +
 *          get(std::integral_constant<std::size_t, 1>{}) + get(std::integral_constant<std::size_t, 2>{});
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_METRICS_PROGRESS_OBSERVER_HPP
#define LAHUTA_PIPELINE_METRICS_PROGRESS_OBSERVER_HPP

#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <mutex>
#include <optional>
#include <string>
#include <utility>

#include <spdlog/fmt/fmt.h>
#include <spdlog/logger.h>
#include <spdlog/sinks/sink.h>

#include "pipeline/metrics/run_observer.hpp"

namespace lahuta::pipeline {

struct ProgressObserverConfig {
  std::chrono::milliseconds interval{50};
  std::string label{"progress"};
  std::optional<std::size_t> total_items;
  std::string marker{"progress"};
  std::string color_on;
  std::string color_off;
  std::shared_ptr<spdlog::logger> progress_logger;
};

class ProgressRunObserver final : public IRunObserver {
  friend class ProgressAwareSink;

public:
  explicit ProgressRunObserver(ProgressObserverConfig config = {})
      : output_(std::cerr), progress_logger_(std::move(config.progress_logger)),
        marker_(std::move(config.marker)),
        label_(config.label.empty() ? "progress" : std::move(config.label)),
        color_on_(std::move(config.color_on)), color_off_(std::move(config.color_off)),
        total_items_(config.total_items), enabled_(config.interval.count() > 0) {
    auto interval = config.interval;
    if (interval < std::chrono::milliseconds{0}) interval = std::chrono::milliseconds{0};
    interval_ns_   = std::chrono::duration_cast<std::chrono::nanoseconds>(interval).count();
    const auto now = Clock::now();
    next_report_ns_.store(to_ns(now) + interval_ns_, std::memory_order_relaxed);
    colorize_ = !color_on_.empty() && !color_off_.empty();
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
      start_time_              = now;
      last_rate_sample_time_   = now;
      last_rate_sample_count_  = 0;
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
    auto lock        = lock_output();
    const auto count = processed_.load(std::memory_order_relaxed);
    if (started_.load(std::memory_order_relaxed) || count > 0) {
      report_locked(count, Clock::now());
      write_raw("\n");
    } else {
      clear_line_locked(true);
    }
    last_len_ = 0;
    last_line_.clear();
  }

  std::uint64_t processed() const noexcept { return processed_.load(std::memory_order_relaxed); }

private:
  using Clock                               = std::chrono::steady_clock;
  static constexpr double RateSampleSeconds = 0.5;

  static std::int64_t to_ns(Clock::time_point tp) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(tp.time_since_epoch()).count();
  }

  void maybe_report(std::uint64_t count) {
    if (interval_ns_ <= 0) return;
    const auto now    = Clock::now();
    const auto now_ns = to_ns(now);
    auto next_ns      = next_report_ns_.load(std::memory_order_relaxed);
    if (now_ns < next_ns) return;
    if (!next_report_ns_.compare_exchange_strong(next_ns, now_ns + interval_ns_, std::memory_order_relaxed))
      return;
    report(count);
  }

  void report(std::uint64_t count) {
    std::lock_guard<std::mutex> lock(output_mutex_);
    report_locked(count, Clock::now());
  }

  void report_locked(std::uint64_t count, Clock::time_point now) {
    const auto skipped    = skipped_.load(std::memory_order_relaxed);
    const auto ok         = count >= skipped ? count - skipped : 0;
    const auto inflight   = inflight_.load(std::memory_order_relaxed);
    const double elapsed  = std::chrono::duration<double>(now - start_time_).count();
    const double avg_rate = elapsed > 0.0 ? static_cast<double>(count) / elapsed : 0.0;

    const auto rate_dt = std::chrono::duration<double>(now - last_rate_sample_time_).count();
    if (rate_dt >= RateSampleSeconds) {
      if (count >= last_rate_sample_count_) {
        last_rate_items_per_sec_ = static_cast<double>(count - last_rate_sample_count_) / rate_dt;
      }
      last_rate_sample_count_ = count;
      last_rate_sample_time_  = now;
    }
    const double rate = last_rate_items_per_sec_ > 0.0 ? last_rate_items_per_sec_ : avg_rate;

    buffer_.clear();
    if (!marker_.empty()) {
      append_colored(buffer_, "[{}]", marker_);
      fmt::format_to(std::back_inserter(buffer_), " ");
    }
    fmt::format_to(std::back_inserter(buffer_), "{}: done=", label_);
    append_colored(buffer_, "{}", count);
    fmt::format_to(std::back_inserter(buffer_), " ok=");
    append_colored(buffer_, "{}", ok);
    fmt::format_to(std::back_inserter(buffer_), " skip=");
    append_colored(buffer_, "{}", skipped);
    fmt::format_to(std::back_inserter(buffer_), " inflight=");
    append_colored(buffer_, "{}", inflight);
    fmt::format_to(std::back_inserter(buffer_), " rate=");
    append_colored(buffer_, "{:.2f}/s", rate);
    fmt::format_to(std::back_inserter(buffer_), " elapsed=");
    append_color_on(buffer_);
    append_duration(buffer_, elapsed);
    append_color_off(buffer_);

    if (total_items_ && *total_items_ > 0) {
      const auto total = *total_items_;
      const double pct = static_cast<double>(count) * 100.0 / static_cast<double>(total);
      fmt::format_to(std::back_inserter(buffer_), " total=");
      append_colored(buffer_, "{}", total);
      fmt::format_to(std::back_inserter(buffer_), " (");
      append_colored(buffer_, "{:.1f}%", pct);
      fmt::format_to(std::back_inserter(buffer_), ")");
      if (rate > 0.0 && count < total) {
        const double eta = static_cast<double>(total - count) / rate;
        fmt::format_to(std::back_inserter(buffer_), " eta=");
        append_color_on(buffer_);
        append_duration(buffer_, eta);
        append_color_off(buffer_);
      }
    }

    last_line_.assign(buffer_.data(), buffer_.size());
    write_line_.clear();
    write_line_.push_back('\r');
    write_line_.append(last_line_);
    write_raw(write_line_);
    last_len_ = last_line_.size();
  }

  void write_raw(std::string_view text) {
    if (progress_logger_) {
      progress_logger_->log(spdlog::level::info, spdlog::string_view_t(text.data(), text.size()));
      progress_logger_->flush();
      return;
    }
    output_ << text << std::flush;
  }

  static void append_duration(fmt::memory_buffer &buf, double seconds) {
    if (seconds < 0.0) seconds = 0.0;
    const auto total_ms = static_cast<std::uint64_t>(seconds * 1000.0);
    const auto ms       = total_ms % 1000;
    const auto total_s  = total_ms / 1000;
    const auto s        = total_s % 60;
    const auto m        = (total_s / 60) % 60;
    const auto h        = total_s / 3600;
    if (h > 0) {
      fmt::format_to(std::back_inserter(buf), "{:02}:{:02}:{:02}", h, m, s);
    } else {
      fmt::format_to(std::back_inserter(buf), "{:02}:{:02}.{:01}", m, s, ms / 100);
    }
  }

  template <typename... Args>
  void append_colored(fmt::memory_buffer &buf, fmt::format_string<Args...> format, Args &&...args) const {
    if (colorize_) {
      buf.append(color_on_);
    }
    fmt::format_to(std::back_inserter(buf), format, std::forward<Args>(args)...);
    if (colorize_) {
      buf.append(color_off_);
    }
  }

  void append_color_on(fmt::memory_buffer &buf) const {
    if (colorize_) {
      buf.append(color_on_);
    }
  }

  void append_color_off(fmt::memory_buffer &buf) const {
    if (colorize_) {
      buf.append(color_off_);
    }
  }

  void clear_line() {
    auto lock = lock_output();
    clear_line_locked(true);
  }

  std::unique_lock<std::mutex> lock_output() { return std::unique_lock<std::mutex>(output_mutex_); }

  void clear_line_locked(bool reset) {
    if (last_len_ == 0) return;
    write_line_.clear();
    write_line_.push_back('\r');
    write_line_.append(last_len_, ' ');
    write_line_.push_back('\r');
    write_raw(write_line_);
    if (reset) {
      last_len_ = 0;
      last_line_.clear();
    }
  }

  void reprint_locked() {
    if (last_len_ == 0 || last_line_.empty()) return;
    write_line_.clear();
    write_line_.push_back('\r');
    write_line_.append(last_line_);
    write_raw(write_line_);
  }

  std::atomic<std::uint64_t> processed_{0};
  std::atomic<std::uint64_t> skipped_{0};
  std::atomic<std::int64_t> inflight_{0};
  std::atomic<std::int64_t> next_report_ns_{0};
  std::int64_t interval_ns_{0};
  std::ostream &output_;
  std::shared_ptr<spdlog::logger> progress_logger_;
  std::string marker_;
  std::string label_;
  std::string color_on_;
  std::string color_off_;
  std::optional<std::size_t> total_items_;
  const bool enabled_{true};
  bool colorize_{false};
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
  std::string write_line_;
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

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_METRICS_PROGRESS_OBSERVER_HPP
