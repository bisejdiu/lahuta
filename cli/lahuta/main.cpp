#include "cli/global_flags.hpp"
#include "cli/global_options.hpp"
#include "logging.hpp"

using namespace lahuta::cli;

// clang-format off
int main(int argc, char* argv[]) {
  try {
    // get global flags and clean argv for subcommands
    auto g = extract_global_flags(argc, argv);

    lahuta::Logger::get_instance().set_log_level(g.log_level);

    if (g.help_requested && g.tail.empty()) {
      print_global_help();
      return 0;
    }

    if (g.tail.empty()) {
      lahuta::Logger::get_logger()->error("No subcommand provided (run lahuta --help for usage)");
      return 1;
    }

    if (g.help_requested) {
      print_global_help();
      return 0;
    }

    const std::string_view subcommand{g.tail[0]};

    if (subcommand.empty() || subcommand.front() == '-') {
      lahuta::Logger::get_logger()->error("Expected subcommand before '{}'. Run lahuta --help for usage.", subcommand);
      return 1;
    }

    const auto& registry = get_command_registry();
    const auto it = registry.find(std::string{subcommand});

    if (it == registry.end()) {
      lahuta::Logger::get_logger()->error("Unknown subcommand '{}' (run lahuta -h for more information)", subcommand);
      return 1;
    }

    // subcommand args
    auto command = it->second();
    return command->run(static_cast<int>(g.tail.size() - 1), g.tail.data() + 1);

  } catch (const std::exception& e) {
    lahuta::Logger::get_logger()->error("{} (run lahuta -h for more information)", e.what());
    return 1;
  }
}
