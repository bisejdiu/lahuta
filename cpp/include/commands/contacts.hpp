#ifndef LAHUTA_CLI_CONTACTS_HPP
#define LAHUTA_CLI_CONTACTS_HPP

#include "commands/command.hpp"
#include <memory>

namespace lahuta::cli {

namespace contacts_opts {
enum ContactsOptionIndex : unsigned {
  Unknown,
  Help,

  // Source options
  SourceDirectory,
  SourceVector,
  SourceFileList,
  SourceDatabase,
  Extension,
  Recursive,

  // Compute options
  Provider,
  InteractionType,

  // Output options
  OutputJson,
  OutputText,
  OutputLog,
  NoCompress,

  // Runtime options
  Threads,
  BatchSize
};
} // namespace contacts_opts

class ContactsCommand final : public CliCommand {
public:
  [[nodiscard]] static std::unique_ptr<CliCommand> create();
  int run(int argc, char* argv[]) override;

private:
  ContactsCommand() = default;
};

} // namespace lahuta::cli

#endif // LAHUTA_CLI_CONTACTS_HPP 
