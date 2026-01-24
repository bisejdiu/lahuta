#ifndef LAHUTA_CLI_GLOBAL_FLAGS_HPP
#define LAHUTA_CLI_GLOBAL_FLAGS_HPP

#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "logging/logging.hpp"

// clang-format off
namespace lahuta::cli {

struct GlobalFlags {
  // Default log level is Warn unless overridden by -v
  lahuta::Logger::LogLevel log_level = lahuta::Logger::LogLevel::Info;
  std::size_t progress_ms = 50;
  std::vector<char*> tail;                                             // argv sans globals
  bool help_requested = false;                                         // -h/--help flag
};

inline GlobalFlags& global_flags_storage() {
  static GlobalFlags flags;
  return flags;
}

inline const GlobalFlags& get_global_flags() {
  return global_flags_storage();
}

inline void set_global_flags(GlobalFlags flags) {
  global_flags_storage() = std::move(flags);
}

inline GlobalFlags extract_global_flags(int argc, char* argv[]) {
  GlobalFlags g;
  g.tail.reserve(argc);

  auto parse_progress_ms = [](std::string_view value) -> std::size_t {
    if (value.empty()) throw std::runtime_error("Option '--progress-ms' expects <ms> (0 disables)");
    std::size_t ms = 0;
    try {
      ms = std::stoull(std::string(value));
    } catch (const std::exception&) {
      throw std::runtime_error("Invalid --progress-ms value '" + std::string(value) + "'");
    }
    return ms;
  };

  int first_non_global = argc;
  for (int i = 1; i < argc; ++i) {
    std::string_view arg{argv[i]};
    if (arg == "-v" || arg == "--verbose") {
      if (i + 1 >= argc) throw std::runtime_error("Option '-v/--verbose' expects <level> (0,1,2)");
      ++i; // Skip level token for now, actual parsing occurs below
      continue;
    }
    if (arg == "--progress-ms") {
      if (i + 1 >= argc) throw std::runtime_error("Option '--progress-ms' expects <ms> (0 disables)");
      ++i;
      continue;
    }
    if (arg.rfind("--progress-ms=", 0) == 0) {
      continue;
    }
    if (arg == "-h" || arg == "--help") {
      g.help_requested = true;
      continue;
    }
    first_non_global = i;
    break;
  }

  for (int i = 1; i < argc; ++i) {
    std::string_view arg{argv[i]};
    if (arg == "-v" || arg == "--verbose") {
      if (i + 1 >= argc) throw std::runtime_error("Option '-v/--verbose' expects <level> (0,1,2)");
      std::string_view lvl{argv[++i]};
      if      (lvl == "0") g.log_level = lahuta::Logger::LogLevel::Error;
      else if (lvl == "1") g.log_level = lahuta::Logger::LogLevel::Info;
      else if (lvl == "2") g.log_level = lahuta::Logger::LogLevel::Debug;
      else throw std::runtime_error("Invalid verbosity level '" + std::string(lvl) + "'. Must be 0, 1 or 2");
      continue;
    } else if (arg == "--progress-ms") {
      if (i + 1 >= argc) throw std::runtime_error("Option '--progress-ms' expects <ms> (0 disables)");
      g.progress_ms = parse_progress_ms(argv[++i]);
      continue;
    } else if (arg.rfind("--progress-ms=", 0) == 0) {
      g.progress_ms = parse_progress_ms(arg.substr(14));
      continue;
    } else if (arg == "-h" || arg == "--help") {
      if (i < first_non_global) continue;
    }
    if (i >= first_non_global) g.tail.push_back(argv[i]);
  }

  set_global_flags(g);
  return g;
}

} // namespace lahuta::cli

#endif // LAHUTA_CLI_GLOBAL_FLAGS_HPP
