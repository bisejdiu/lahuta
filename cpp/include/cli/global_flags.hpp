#ifndef LAHUTA_CLI_GLOBAL_FLAGS_HPP
#define LAHUTA_CLI_GLOBAL_FLAGS_HPP

#include <vector>
#include <string_view>
#include <stdexcept>
#include "logging.hpp"

namespace lahuta::cli {

struct GlobalFlags {
  lahuta::Logger::LogLevel log_level = lahuta::Logger::LogLevel::Warn; // default 1
  std::vector<char*> tail;                                             // argv sans globals
  bool help_requested = false;                                         // -h/--help flag
};

inline GlobalFlags extract_global_flags(int argc, char* argv[]) {
  GlobalFlags g;
  g.tail.reserve(argc);

  // check subcommand
  bool has_explicit_subcommand = false;
  for (int i = 1; i < argc; ++i) {
    std::string_view arg{argv[i]};
    if (!arg.empty() && arg[0] != '-') {
      if (arg == "contacts") { has_explicit_subcommand = true; }
      break;
    }
  }

  for (int i = 1; i < argc; ++i) { // skip executable name
    std::string_view arg{argv[i]};
    if (arg == "-v" || arg == "--verbose") {
      if (i + 1 >= argc) throw std::runtime_error("Option '-v/--verbose' expects <level> (0,1,2)");

      std::string_view lvl{argv[++i]};
      if      (lvl == "0") g.log_level = lahuta::Logger::LogLevel::Error;
      else if (lvl == "1") g.log_level = lahuta::Logger::LogLevel::Warn;
      else if (lvl == "2") g.log_level = lahuta::Logger::LogLevel::Debug;
      else throw std::runtime_error("Invalid verbosity level '" + std::string(lvl) + "'. Must be 0, 1 or 2");

      continue;           // swallow both tokens
    } else if (arg == "-h" || arg == "--help") {
      // only treat help as global if there's no explicit subcommand
      if (!has_explicit_subcommand) {
        g.help_requested = true;
        continue;           // swallow help flag
      }
      // let the subcommand handle it
    }
    g.tail.push_back(argv[i]);
  }

  // if no subcommand was explicitly provided, default to "contacts"
  if (!has_explicit_subcommand && !g.tail.empty()) {
    g.tail.insert(g.tail.begin(), const_cast<char*>("contacts"));
  } else if (!has_explicit_subcommand && g.tail.empty()) {
    g.tail.push_back(const_cast<char*>("contacts"));
  }

  return g;
}

} // namespace lahuta::cli

#endif // LAHUTA_CLI_GLOBAL_FLAGS_HPP
