#ifndef LAHUTA_CLI_GLOBAL_ARGS_HPP
#define LAHUTA_CLI_GLOBAL_ARGS_HPP

#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "specs/command_spec.hpp"

namespace lahuta::cli {

struct GlobalArgSplit {
  std::vector<const char *> global_args;
  std::vector<const char *> tail;
};

inline std::string build_global_usage(const CommandSpecRegistry &registry) {
  std::size_t max_name_len = 0;
  for (const auto &[name, spec] : registry) {
    if (spec) max_name_len = std::max(max_name_len, name.size());
  }

  std::string usage = "Usage: lahuta [global-options] <subcommand> [subcommand-options]\n\n"
                      "Lahuta runs structural analysis, such as contacts computations, at scale.\n\n"
                      "Available subcommands:\n";

  const std::size_t name_pad = max_name_len + 2;
  for (const auto &[name, spec] : registry) {
    if (!spec) continue;

    usage.append("  ");
    usage.append(name);
    if (name_pad > name.size()) {
      usage.append(name_pad - name.size(), ' ');
    }
    usage.append(spec->summary());
    usage.push_back('\n');
  }

  usage.append("\nGlobal Options:");
  return usage;
}

inline GlobalArgSplit split_global_args(int argc, char *argv[]) {
  GlobalArgSplit split;
  if (argc <= 0) return split;

  split.global_args.reserve(static_cast<std::size_t>(argc));
  split.tail.reserve(static_cast<std::size_t>(argc));
  split.global_args.push_back(argv[0]);

  int first_non_global = argc;
  for (int i = 1; i < argc; ++i) {
    const std::string_view arg{argv[i]};
    if (arg == "-v" || arg == "--verbose") {
      if (i + 1 >= argc) {
        throw std::runtime_error("Option '-v/--verbose' expects <level> (0,1,2)");
      }
      ++i;
      continue;
    }
    if (arg == "--progress-ms") {
      if (i + 1 >= argc) {
        throw std::runtime_error("Option '--progress-ms' expects <ms> (0 disables)");
      }
      ++i;
      continue;
    }
    if (arg.rfind("--progress-ms=", 0) == 0) {
      continue;
    }
    if (arg == "--progress-no-color") {
      continue;
    }
    if (arg == "-h" || arg == "--help") {
      continue;
    }
    if (arg == "--version") {
      continue;
    }
    first_non_global = i;
    break;
  }

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg{argv[i]};
    if (arg == "-v" || arg == "--verbose") {
      if (i + 1 >= argc) {
        throw std::runtime_error("Option '-v/--verbose' expects <level> (0,1,2)");
      }
      split.global_args.push_back(argv[i]);
      split.global_args.push_back(argv[i + 1]);
      ++i;
      continue;
    }
    if (arg == "--progress-ms") {
      if (i + 1 >= argc) {
        throw std::runtime_error("Option '--progress-ms' expects <ms> (0 disables)");
      }
      split.global_args.push_back(argv[i]);
      split.global_args.push_back(argv[i + 1]);
      ++i;
      continue;
    }
    if (arg.rfind("--progress-ms=", 0) == 0) {
      split.global_args.push_back(argv[i]);
      continue;
    }
    if (arg == "--progress-no-color") {
      split.global_args.push_back(argv[i]);
      continue;
    }
    if (arg == "-h" || arg == "--help") {
      if (i < first_non_global) {
        split.global_args.push_back(argv[i]);
        continue;
      }
    }
    if (arg == "--version") {
      if (i < first_non_global) {
        split.global_args.push_back(argv[i]);
        continue;
      }
    }

    if (i >= first_non_global) {
      split.tail.push_back(argv[i]);
    }
  }

  return split;
}

} // namespace lahuta::cli

#endif // LAHUTA_CLI_GLOBAL_ARGS_HPP
