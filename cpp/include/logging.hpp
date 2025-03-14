#ifndef LAHUTA_LOGGING_HPP
#define LAHUTA_LOGGING_HPP

#include "spinner.hpp"
#include <spdlog/common.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#define LAHUTA_TRACE EntryExitTraceLogger scoped_logger(__FUNCTION__, __FILE__, __LINE__)

// clang-format off
namespace lahuta {

class Logger {
public:
  enum class LogLevel    { Trace, Debug, Info, Warn, Error, Critical, Off };
  enum class FormatStyle { Simple, Detailed };

  static Logger &get_instance() {
    static Logger instance;
    return instance;
  }

  void set_format(FormatStyle style) {
    switch (style) {
      case FormatStyle::Detailed:
        spdlog::set_pattern("[%T] [thread %t] [%^%l%$] %t] %v");
        break;
      case FormatStyle::Simple:
        spdlog::set_pattern("[%^%l%$] %v");
        break;
    }
  }

  void set_log_level(LogLevel level) { spdlog::set_level(convert_log_level(level)); }
  void log(LogLevel level, const std::string &message) { logger->log(convert_log_level(level), message); }
  LogLevel get_log_level() { return static_cast<LogLevel>(logger->level()); }

    void configure_for_spinner(indicators::MinimalProgressSpinner* spinner, std::mutex &spinner_mutex) {
        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        // wrap the console sink with the spinner-aware sink.
        auto spinner_sink = std::make_shared<indicators::SpinnerAwareSink>(console_sink, spinner, spinner_mutex);
        // create a logger with the spinner-aware sink.
        logger = std::make_shared<spdlog::logger>("console", spinner_sink);
        spdlog::set_default_logger(logger);
    }

private:
  std::shared_ptr<spdlog::logger> logger;

  Logger() {
    logger = spdlog::stdout_color_mt("console");
    set_log_level(LogLevel::Info);
    set_format(FormatStyle::Simple);
  }

  spdlog::level::level_enum convert_log_level(LogLevel level) {
    switch (level) {
      case LogLevel::Trace:    return spdlog::level::trace;
      case LogLevel::Debug:    return spdlog::level::debug;
      case LogLevel::Info:     return spdlog::level::info;
      case LogLevel::Warn:     return spdlog::level::warn;
      case LogLevel::Error:    return spdlog::level::err;
      case LogLevel::Critical: return spdlog::level::critical;
      case LogLevel::Off:      return spdlog::level::off;
    }
    return spdlog::level::info;
  }
};


class EntryExitTraceLogger {
public:
  EntryExitTraceLogger(const std::string &func, const char *file, int line)
      : func_name(func), file(file), line(line) {
    spdlog::trace("Entering {} ({}:{})", func_name, file, line);
  }
  ~EntryExitTraceLogger() { spdlog::trace("Exiting {} ({}:{})", func_name, file, line); }

private:
  std::string func_name;
  const char *file;
  int line;
};

} // namespace lahuta

#endif // LAHUTA_LOGGING_HPP
