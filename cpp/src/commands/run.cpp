#include "commands/run.hpp"
#include "cli/arg_validation.hpp"
#include "contacts/arpeggio/provider.hpp"
#include "contacts/engine.hpp"
#include "contacts/molstar/provider.hpp"
#include "entities/formatter.hpp"
#include "lahuta.hpp"
#include "logging.hpp"
#include <iostream>
#include <vector>

// clang-format off
namespace lahuta::cli {

namespace run_opts {
const option::Descriptor usage[] = {
  {RunOptionIndex::Unknown, 0, "", "", validate::Unknown,
   "Usage: lahuta run <input_file> [options]\n\n"
   "Compute inter-atomic contacts from molecular structure files.\n\n"
   "Options:"},
  {RunOptionIndex::Help, 0, "h", "help", option::Arg::None,
   "  --help, -h                   \tPrint this help message and exit."},
  {RunOptionIndex::Provider, 0, "p", "provider", validate::Provider,
   "  --provider, -p <provider>    \tContact provider: 'arpeggio' or 'molstar' (default: molstar)."},
  {RunOptionIndex::ContactType, 0, "t", "type", validate::ContactType,
   "  --type, -t <type>            \tContact type to compute. Can be specified multiple times.\n"
   "                               \tIf not specified, all available types are computed."},
  {RunOptionIndex::Quiet, 0, "q", "quiet", option::Arg::None,
   "  --quiet, -q                  \tSuppress verbose output."}, {0, 0, "", "", option::Arg::None,
   "\nAvailable Contact Types:\n"
   "  Common (both providers):     hbond, hydrophobic, ionic\n"
   "  MolStar specific:            weak_hbond, halogen, metalic, cationpi, pistacking\n"
   "  Arpeggio specific:           polar_hbond, weak_polar_hbond, aromatic, carbonyl,\n"
   "                               vdw, donor_pi, sulphur_pi, carbon_pi"},
  {0, 0, 0, 0, 0, 0}
};
} // namespace run_opts

[[nodiscard]] std::unique_ptr<CliCommand> RunCommand::create() {
    return std::unique_ptr<CliCommand>(new RunCommand());
}

template<typename Provider>
void RunCommand::compute_contacts(const Topology& topology, const std::vector<std::string>& contact_types, bool quiet) {
  InteractionEngine<Provider> engine;

  if (contact_types.empty()) {
    // Compute all available contact types
    auto all_contacts = engine.compute(topology);
    if (!quiet) Logger::get_logger()->info("Computing all available contact types...");
    ContactTableFormatter::print_contact_table(all_contacts, topology, "All Contacts");
  } else {
    // Compute specific contact types
    for (const auto& type_str : contact_types) {
      const auto interaction_type = get_interaction_type(type_str);
      if (interaction_type == InteractionType::None) {
        Logger::get_logger()->warn("Unknown contact type '{}', skipping.", type_str);
        continue;
      }

      if (!quiet) Logger::get_logger()->info("Computing {} contacts...", type_str);

      auto contacts = engine.compute(topology, interaction_type);
      ContactTableFormatter::print_contact_table(contacts, topology, type_str);
    }
  }
}

int RunCommand::run(int argc, char* argv[]) {
  option::Stats stats(true, run_opts::usage, argc, const_cast<const char**>(argv));
  std::vector<option::Option> options(stats.options_max);
  std::vector<option::Option> buffer (stats.buffer_max);
  option::Parser parse(true, run_opts::usage, argc, const_cast<const char**>(argv), options.data(), buffer.data());

  if (parse.error()) return 1;

  if (options[run_opts::RunOptionIndex::Help]) {
    option::printUsage(std::cout, run_opts::usage);
    return 0;
  }

  if (parse.nonOptionsCount() == 0) {
    Logger::get_logger()->error("Input file is required{}", validate::RUN_HELP_MSG_SUFFIX);
    return 1;
  }

  const std::string input_file{parse.nonOption(0)};

  // Parse provider option
  std::string provider = "molstar";
  if (const auto& provider_opt = options[run_opts::RunOptionIndex::Provider]; provider_opt) {
    provider = provider_opt.arg;
  }

  // Parse contact types
  std::vector<std::string> contact_types;
  for (option::Option* opt = options[run_opts::RunOptionIndex::ContactType]; opt; opt = opt->next()) {
    contact_types.emplace_back(opt->arg);
  }

  // Check quiet flag (note: global quiet is handled by the caller)
  const bool quiet = static_cast<bool>(options[run_opts::RunOptionIndex::Quiet]);

  if (quiet) { Logger::get_instance().set_log_level(Logger::LogLevel::Error); }
  else       { Logger::get_instance().set_log_level(Logger::LogLevel::Info); }

  try {

    Logger::get_logger()->info("Loading file: {}", input_file);

    Luni luni(input_file);

    // Set the appropriate contact provider
    const auto computer_type = (provider == "arpeggio") ? ContactComputerType::Arpeggio : ContactComputerType::Molstar;
    luni.set_atom_typing_method(computer_type);

    Logger::get_logger()->info("Building topology with {} provider...", provider);

    if (!luni.build_topology()) {
      Logger::get_logger()->error("Failed to build topology for file: {}{}", input_file, validate::HELP_MSG_SUFFIX);
      return 1;
    }

    const auto& topology = luni.get_topology();
    Logger::get_logger()->info("Molecule loaded: {} atoms, {} bonds", luni.get_molecule().getNumAtoms(), luni.get_molecule().getNumBonds());

    // Compute contacts
    if (provider == "arpeggio") { compute_contacts<ArpeggioContactProvider>(topology, contact_types, quiet); }
    else                        { compute_contacts<MolStarContactProvider> (topology, contact_types, quiet); }

  } catch (const std::exception& e) {
    Logger::get_logger()->error("{}{}", e.what(), validate::HELP_MSG_SUFFIX);
    return 1;
  }

  return 0;
}

} // namespace lahuta::cli 
