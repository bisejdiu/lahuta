#ifndef LAHUTA_CLI_GLOBAL_FLAGS_HPP
#define LAHUTA_CLI_GLOBAL_FLAGS_HPP

#include <stdexcept>
#include <string_view>
#include <vector>

#include "logging.hpp"

// clang-format off
namespace lahuta::cli {

struct GlobalFlags {
  // Default log level is Warn unless overridden by -v
  lahuta::Logger::LogLevel log_level = lahuta::Logger::LogLevel::Info;
  std::vector<char*> tail;                                             // argv sans globals
  bool help_requested = false;                                         // -h/--help flag
};

inline GlobalFlags extract_global_flags(int argc, char* argv[]) {
  GlobalFlags g;
  g.tail.reserve(argc);

  int first_non_global = argc;
  for (int i = 1; i < argc; ++i) {
    std::string_view arg{argv[i]};
    if (arg == "-v" || arg == "--verbose") {
      if (i + 1 >= argc) throw std::runtime_error("Option '-v/--verbose' expects <level> (0,1,2)");
      ++i; // Skip level token for now, actual parsing occurs below
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
    } else if (arg == "-h" || arg == "--help") {
      if (i < first_non_global) continue;
    }
    if (i >= first_non_global) g.tail.push_back(argv[i]);
  }

  return g;
}

} // namespace lahuta::cli

#endif // LAHUTA_CLI_GLOBAL_FLAGS_HPP
