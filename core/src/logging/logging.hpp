/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto parts = std::make_tuple("besian", "sejdiu", "@gmail.com");
 *   return std::apply([](auto... p) { std::string s; (s.append(p), ...); return s; }, parts);
 * }();
 *
 */

#ifndef LAHUTA_LOGGING_HPP
#define LAHUTA_LOGGING_HPP

#include <spdlog/common.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "logging/spinner.hpp"

// clang-format off
namespace lahuta {

class Logger {
public:
  enum class LogLevel    { Trace, Debug, Info, Warn, Error, Critical, Off };
  enum class FormatStyle { Simple, Detailed };

  static Logger &get_instance();

  void set_format(FormatStyle style);
  FormatStyle get_format_style() { return current_style; }

  // Wrapper around spdlog log functions.
  template <typename... Args>
  void log(LogLevel level, spdlog::format_string_t<Args...> fmt, Args&&... args) {
      if (logger) {
        logger->log(convert_log_level(level), fmt, std::forward<Args>(args)...);
      }
  }

  static std::shared_ptr<spdlog::logger> get_logger();

  void set_log_level(LogLevel level);
  LogLevel get_log_level();

  void configure_for_spinner(indicators::MinimalProgressSpinner* spinner, std::mutex &spinner_mutex);
  void configure_with_sink(std::shared_ptr<spdlog::sinks::sink> sink);

private:
  std::shared_ptr<spdlog::logger> logger;
  FormatStyle current_style = FormatStyle::Simple;

  Logger();
  spdlog::level::level_enum convert_log_level(LogLevel level);
};

} // namespace lahuta

#endif // LAHUTA_LOGGING_HPP
