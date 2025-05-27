#ifndef LAHUTA_CLI_RUN_HPP
#define LAHUTA_CLI_RUN_HPP

#include "commands/command.hpp"
#include <memory>
#include <string>
#include <vector>

// clang-format off
namespace lahuta {
class Topology;

namespace cli {

namespace run_opts {
enum RunOptionIndex : unsigned {
  Unknown,
  Help,
  Provider,
  ContactType,
  Quiet
};
} // namespace run_opts

class RunCommand final : public CliCommand {
public:
  [[nodiscard]] static std::unique_ptr<CliCommand> create();
  int run(int argc, char* argv[]) override;

private:
  RunCommand() = default;

  // computes contacts given a provider and a list of contact types
  template<typename Provider>
  static void compute_contacts(const Topology& topology, const std::vector<std::string>& contact_types, bool quiet);
};

} // namespace cli
} // namespace lahuta

#endif // LAHUTA_CLI_RUN_HPP
