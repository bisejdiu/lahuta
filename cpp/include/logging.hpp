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

  void set_format(FormatStyle style);

  // Wrapper around spdlog log functions.
  template <typename... Args>
  void log(LogLevel level, spdlog::format_string_t<Args...> fmt, Args&&... args) {
      spdlog::log(convert_log_level(level), fmt, std::forward<Args>(args)...);
  }

  static std::shared_ptr<spdlog::logger> get_logger();

  void set_log_level(LogLevel level);
  LogLevel get_log_level();

  void configure_for_spinner(indicators::MinimalProgressSpinner* spinner, std::mutex &spinner_mutex);

private:
  std::shared_ptr<spdlog::logger> logger;

  Logger();
  spdlog::level::level_enum convert_log_level(LogLevel level);
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
