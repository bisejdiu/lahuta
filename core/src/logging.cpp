#include "logging.hpp"

namespace lahuta {

Logger &Logger::get_instance() {
  static Logger instance;
  return instance;
}

void Logger::set_format(FormatStyle style) {
  current_style = style;
  if (logger) {
    switch (style) {
      case FormatStyle::Detailed:
        logger->set_pattern("[%T] [thread %t] [%^%l%$] %t] %v");
        break;
      case FormatStyle::Simple:
        logger->set_pattern("[%^%l%$] %v");
        break;
    }
  }
}

std::shared_ptr<spdlog::logger> Logger::get_logger() { return get_instance().logger; }

void Logger::set_log_level(LogLevel level) {
  if (logger) logger->set_level(convert_log_level(level));
}

Logger::LogLevel Logger::get_log_level() { return static_cast<LogLevel>(logger->level()); }

void Logger::configure_for_spinner(indicators::MinimalProgressSpinner *spinner, std::mutex &spinner_mutex) {
  auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  std::shared_ptr<spdlog::sinks::sink> sink;
  if (spinner) {
    sink = std::make_shared<indicators::SpinnerAwareSink>(console_sink, spinner, spinner_mutex);
  } else {
    sink = console_sink;
  }
  auto new_logger   = std::make_shared<spdlog::logger>("console", sink);

  new_logger->set_level(logger->level());
  logger = new_logger;
  spdlog::set_default_logger(new_logger);
}

Logger::Logger() {
  auto existing = spdlog::get("console");
  if (existing) {
    logger = existing;
  } else {
    logger = spdlog::stdout_color_mt("console");
  }
  // Make this the default logger within this registry to keep free-function calls consistent in this TU
  if (logger) {
    spdlog::set_default_logger(logger);
  }
  set_log_level(LogLevel::Info);
  set_format(FormatStyle::Simple);
}

spdlog::level::level_enum Logger::convert_log_level(LogLevel level) {
    switch (level) {
        case LogLevel::Trace:    return spdlog::level::trace;
        case LogLevel::Debug:    return spdlog::level::debug;
        case LogLevel::Info:     return spdlog::level::info;
        case LogLevel::Warn:     return spdlog::level::warn;
        case LogLevel::Error:    return spdlog::level::err;
        case LogLevel::Critical: return spdlog::level::critical;
        case LogLevel::Off:      return spdlog::level::off;
    }
}

} // namespace lahuta
