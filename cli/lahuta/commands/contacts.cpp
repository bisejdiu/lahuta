#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <variant>
#include <vector>

#include "analysis/contacts/computation.hpp"
#include "analysis/contacts/provider.hpp"
#include "cli/arg_validation.hpp"
#include "cli/extension_utils.hpp"
#include "commands/contacts.hpp"
#include "db/db.hpp"
#include "gemmi/third_party/stb_sprintf.h"
#include "io/sinks/logging.hpp"
#include "io/sinks/ndjson.hpp"
#include "logging.hpp"
#include "pipeline/compute/parameters.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/dynamic/sources.hpp"
#include "runtime.hpp"

// clang-format off
namespace lahuta::cli {

using namespace lahuta::pipeline;

using Source = std::variant<sources::Directory, std::vector<std::string>, sources::FileList>;

static void initialize_runtime(int num_threads) {
  lahuta::LahutaRuntime::ensure_initialized(static_cast<std::size_t>(num_threads));
}

struct ContactsOptions {
  enum class SourceMode { Directory, Vector, FileList, Database };

  SourceMode source_mode = SourceMode::Directory;
  std::string directory_path;
  std::vector<std::string> extensions{".cif", ".cif.gz"};
  bool recursive = true;
  std::vector<std::string> file_vector;
  std::string file_list_path;
  std::string database_path;

  analysis::contacts::ContactProvider provider = analysis::contacts::ContactProvider::MolStar;
  InteractionTypeSet interaction_types = InteractionTypeSet::all();

  bool want_json   = false;
  bool want_text   = false;
  bool want_log    = false;
  bool no_compress = false;

  int threads = 8;
  size_t batch_size = 200;
  size_t writer_threads = 1;
};

Source pick_source(const ContactsOptions& cli) {
  switch (cli.source_mode) {
    case ContactsOptions::SourceMode::Directory: return sources::Directory{cli.directory_path, cli.extensions, cli.recursive, cli.batch_size};
    case ContactsOptions::SourceMode::Vector:    return cli.file_vector;
    case ContactsOptions::SourceMode::FileList:  return sources::FileList {cli.file_list_path};
    case ContactsOptions::SourceMode::Database:  break; // Handled separately
  }
  throw std::logic_error("Invalid source mode");
}

namespace contacts_opts {
const option::Descriptor usage[] = {
  {ContactsOptionIndex::Unknown, 0, "", "", validate::Unknown,
   "Usage: lahuta contacts [options]\n\n"
   "Compute inter-atomic contacts.\n\n"
   "Input Options (choose one):"},
  {ContactsOptionIndex::Help, 0, "h", "help", option::Arg::None,
   "  --help, -h                   \tPrint this help message and exit."},
  {ContactsOptionIndex::SourceDirectory, 0, "d", "directory", validate::Required,
   "  --directory, -d <path>       \tProcess all files in directory."},
  {ContactsOptionIndex::SourceVector, 0, "f", "files", validate::Required,
   "  --files, -f <file1,file2>    \tProcess specific files (comma-separated or repeat -f)."},
  {ContactsOptionIndex::SourceFileList, 0, "l", "file-list", validate::Required,
   "  --file-list, -l <path>       \tProcess files listed in text file (one per line)."},
  {ContactsOptionIndex::SourceDatabase, 0, "", "database", validate::Required,
   "  --database <path>            \tProcess structures from database."},
  {ContactsOptionIndex::Extension, 0, "e", "extension", validate::Required,
   "  --extension, -e <ext>        \tFile extension(s) for directory mode. Repeat or comma-separate values (default: .cif, .cif.gz)."},
  {ContactsOptionIndex::Recursive, 0, "r", "recursive", option::Arg::None,
   "  --recursive, -r              \tRecursively search subdirectories."},
  {0, 0, "", "", option::Arg::None,
   "\nCompute Options:"},
  {ContactsOptionIndex::Provider, 0, "p", "provider", validate::Provider,
   "  --provider, -p <provider>    \tContact provider: 'molstar', 'arpeggio', or 'getcontacts' (default: molstar)."},
  {ContactsOptionIndex::InteractionType, 0, "i", "interaction", validate::ContactType,
   "  --interaction, -i <type>     \tInteraction type(s): 'hbond', 'hydrophobic', 'ionic', etc. Repeat or comma-separate to combine.\n"},
  {0, 0, "", "", option::Arg::None,
   "\nOutput Options:"},
  {ContactsOptionIndex::OutputJson, 0, "", "json", option::Arg::None,
   "  --json                       \tOutput results in JSON format."},
  {ContactsOptionIndex::OutputText, 0, "", "text", option::Arg::None,
   "  --text                       \tOutput results in text format."},
  {ContactsOptionIndex::OutputLog, 0, "", "log", option::Arg::None,
   "  --log                        \tOutput results to standard output (logging)."},
  {ContactsOptionIndex::NoCompress, 0, "", "no-compress", option::Arg::None,
   "  --no-compress                \tDisable gzip compression (default: enabled)."},
  {0, 0, "", "", option::Arg::None,
   "\nRuntime Options:"},
  {ContactsOptionIndex::Threads, 0, "t", "threads", validate::Required,
   "  --threads, -t <num>          \tNumber of threads to use (default: 8)."},
  {ContactsOptionIndex::BatchSize, 0, "b", "batch-size", validate::Required,
   "  --batch-size, -b <size>      \tBatch size for processing (default: 200)."},
  {ContactsOptionIndex::WriterThreads, 0, "", "writer-threads", validate::Required,
   "  --writer-threads <num>       \tNumber of writer threads per sink (default: 1)."},
  {0, 0, 0, 0, 0, 0}
};
} // namespace contacts_opts

[[nodiscard]] std::unique_ptr<CliCommand> ContactsCommand::create() {
  return std::unique_ptr<CliCommand>(new ContactsCommand());
}

int ContactsCommand::run(int argc, char* argv[]) {
  option::Stats stats(true, contacts_opts::usage, argc, const_cast<const char**>(argv));
  std::vector<option::Option> options(stats.options_max);
  std::vector<option::Option> buffer (stats.buffer_max);
  option::Parser parse(true, contacts_opts::usage, argc, const_cast<const char**>(argv), options.data(), buffer.data());

  if (parse.error()) return 1;

  if (options[contacts_opts::ContactsOptionIndex::Help]) {
    option::printUsage(std::cout, contacts_opts::usage);
    return 0;
  }

  try {
    // parse CLI Arguments into ContactsOptions
    ContactsOptions cli;

    // parse source options
    int source_count = 0;
    if (options[contacts_opts::ContactsOptionIndex::SourceDirectory]) {
      cli.source_mode = ContactsOptions::SourceMode::Directory;
      cli.directory_path = options[contacts_opts::ContactsOptionIndex::SourceDirectory].arg;
      source_count++;
    }
    if (options[contacts_opts::ContactsOptionIndex::SourceVector]) {
      cli.source_mode = ContactsOptions::SourceMode::Vector;
      cli.file_vector.clear();
      for (const option::Option* opt = &options[contacts_opts::ContactsOptionIndex::SourceVector];
           opt != nullptr;
           opt = opt->next()) {
        if (opt->arg) parse_file_argument(opt->arg, cli.file_vector);
      }
      source_count++;
    }
    if (options[contacts_opts::ContactsOptionIndex::SourceFileList]) {
      cli.source_mode = ContactsOptions::SourceMode::FileList;
      cli.file_list_path = options[contacts_opts::ContactsOptionIndex::SourceFileList].arg;
      source_count++;
    }
    if (options[contacts_opts::ContactsOptionIndex::SourceDatabase]) {
      cli.source_mode = ContactsOptions::SourceMode::Database;
      cli.database_path = options[contacts_opts::ContactsOptionIndex::SourceDatabase].arg;
      source_count++;
    }

    if (source_count == 0) {
      Logger::get_logger()->error("Must specify exactly one source option: --directory, --files, --file-list, or --database");
      return 1;
    }
    if (source_count > 1) {
      Logger::get_logger()->error("Cannot specify multiple source options");
      return 1;
    }

    // Parse other options
    if (options[contacts_opts::ContactsOptionIndex::Extension]) {
      cli.extensions.clear();
      for (const option::Option* opt = &options[contacts_opts::ContactsOptionIndex::Extension];
           opt != nullptr;
           opt = opt->next()) {
        parse_extension_argument(opt->arg ? opt->arg : "", cli.extensions);
      }
    }
    if (cli.extensions.empty()) cli.extensions.emplace_back();

    cli.recursive = options[contacts_opts::ContactsOptionIndex::Recursive] ? true : false;

    // Parse provider
    if (options[contacts_opts::ContactsOptionIndex::Provider]) {
      std::string_view provider = options[contacts_opts::ContactsOptionIndex::Provider].arg;
      if (auto p = analysis::contacts::contact_provider_from_string(provider)) {
        cli.provider = *p;
      } else {
        Logger::get_logger()->error("Invalid provider '{}'. Must be 'molstar', 'arpeggio', or 'getcontacts'", provider);
        return 1;
      }
    }

    // Parse interaction type
    using coi = contacts_opts::ContactsOptionIndex;
    if (options[coi::InteractionType]) {
      bool has_selection = false;
      InteractionTypeSet selected;
      for (const option::Option* opt = &options[coi::InteractionType]; opt != nullptr; opt = opt->next()) {
        std::string_view arg = opt->arg ? opt->arg : "";
        auto parsed = parse_interaction_type_sequence(arg, ',');
        if (!parsed) parsed = parse_interaction_type_sequence(arg, '|');
        if (!parsed) {
          Logger::get_logger()->error("Invalid interaction type specification '{}'", arg);
          return 1;
        }
        if (!has_selection) {
          selected = *parsed;
          has_selection = true;
        } else {
          selected |= *parsed;
        }
      }
      if (has_selection) cli.interaction_types = selected;
    }

    // Parse output options
    cli.want_json   = options[contacts_opts::ContactsOptionIndex::OutputJson] ? true : false;
    cli.want_text   = options[contacts_opts::ContactsOptionIndex::OutputText] ? true : false;
    cli.want_log    = options[contacts_opts::ContactsOptionIndex::OutputLog]  ? true : false;
    cli.no_compress = options[contacts_opts::ContactsOptionIndex::NoCompress] ? true : false;

    if (!cli.want_json && !cli.want_text && !cli.want_log) cli.want_json = true; // json if nothing specified

    // runtime options
    if (options[contacts_opts::ContactsOptionIndex::Threads]) {
      cli.threads = std::stoi(options[contacts_opts::ContactsOptionIndex::Threads].arg);
      if (cli.threads <= 0) {
        Logger::get_logger()->error("Threads must be positive");
        return 1;
      }
    }

    if (options[contacts_opts::ContactsOptionIndex::BatchSize]) {
      cli.batch_size = std::stoull(options[contacts_opts::ContactsOptionIndex::BatchSize].arg);
      if (cli.batch_size == 0) {
        Logger::get_logger()->error("Batch size must be positive");
        return 1;
      }
    }

    if (options[contacts_opts::ContactsOptionIndex::WriterThreads]) {
      cli.writer_threads = std::stoull(options[contacts_opts::ContactsOptionIndex::WriterThreads].arg);
      if (cli.writer_threads == 0) {
        Logger::get_logger()->error("Writer threads must be positive");
        return 1;
      }
    }

    initialize_runtime(cli.threads);

    if (cli.source_mode == ContactsOptions::SourceMode::Directory) {
      Logger::get_logger()->info("Source directory: {}", cli.directory_path);
      Logger::get_logger()->info("Extensions: {}", describe_extensions(cli.extensions));
      Logger::get_logger()->info("Recursive: {}", cli.recursive ? "Yes" : "No");
    }

    bool is_db = (cli.source_mode == ContactsOptions::SourceMode::Database);
    if (is_db) {
      auto db = std::make_shared<LMDBDatabase>(cli.database_path);
      auto src = dynamic::sources_factory::from_lmdb(db, std::string{}, cli.batch_size);
      dynamic::StageManager mgr(std::move(src));
      // Enable conditional built-ins injection for CLI ergonomics
      mgr.set_auto_builtins(true);

      // Configure built-ins
      mgr.get_system_params().is_model = true;
      mgr.get_topology_params().atom_typing_method = analysis::contacts::typing_for_provider(cli.provider);

      const bool json_out = (cli.want_json || !cli.want_text);
      // Add contacts as a compute-backed task using ContactsKernel
      {
        pipeline::compute::ContactsParams p{};
        p.provider = cli.provider;
        p.type     = cli.interaction_types;
        p.channel  = "contacts";
        p.format   = json_out ? pipeline::compute::ContactsOutputFormat::Json
                              : pipeline::compute::ContactsOutputFormat::Text;
        mgr.add_computation(
          "contacts",
          {},
          [label = std::string("contacts"), p]() {
            return std::make_unique<analysis::contacts::ContactsComputation>(label, p);
          },
          /*thread_safe=*/true
        );
      }

      // Sinks
      dynamic::BackpressureConfig sink_cfg;
      sink_cfg.writer_threads = cli.writer_threads;
      if  (json_out && cli.want_json) mgr.connect_sink("contacts", std::make_shared<dynamic::NdjsonFileSink>("contacts.jsonl"), sink_cfg);
      if (!json_out && cli.want_text) mgr.connect_sink("contacts", std::make_shared<dynamic::NdjsonFileSink>("contacts.txt"), sink_cfg);
      if (cli.want_log)  mgr.connect_sink("contacts", std::make_shared<dynamic::LoggingSink>(), sink_cfg);

      mgr.compile();
      mgr.run(static_cast<std::size_t>(cli.threads));
    } else {
      Source src_variant = pick_source(cli);
      std::visit([&](auto&& src) {
        using SrcT = std::decay_t<decltype(src)>;
        auto source_ptr = std::unique_ptr<lahuta::sources::IDescriptor>{};
        if constexpr (std::is_same_v<SrcT, sources::Directory>) {
          source_ptr = dynamic::sources_factory::from_directory(std::move(src));
        } else if constexpr (std::is_same_v<SrcT, std::vector<std::string>>) {
          source_ptr = dynamic::sources_factory::from_vector(std::move(src));
        } else if constexpr (std::is_same_v<SrcT, sources::FileList>) {
          source_ptr = dynamic::sources_factory::from_filelist(std::move(src));
        } else {
          static_assert(sizeof(SrcT) == 0, "Unsupported source type");
        }
        dynamic::StageManager mgr(std::move(source_ptr));
        // Enable conditional built-ins injection for CLI ergonomics
        mgr.set_auto_builtins(true);

        // Configure built-ins
        mgr.get_system_params().is_model = false;
        mgr.get_topology_params().atom_typing_method = analysis::contacts::typing_for_provider(cli.provider);

        // Tasks: ensure_typing -> contacts
        const bool json_out = (cli.want_json || !cli.want_text);
        {
          pipeline::compute::ContactsParams p{};
          p.provider = cli.provider;
          p.type     = cli.interaction_types;
          p.channel  = "contacts";
          p.format   = json_out ? pipeline::compute::ContactsOutputFormat::Json
                                : pipeline::compute::ContactsOutputFormat::Text;
          mgr.add_computation(
            "contacts",
            {},
            [label = std::string("contacts"), p]() {
              return std::make_unique<analysis::contacts::ContactsComputation>(label, p);
            },
            /*thread_safe=*/true
          );
        }

        // Sinks
        dynamic::BackpressureConfig sink_cfg;
        sink_cfg.writer_threads = cli.writer_threads;
        if  (json_out && cli.want_json) mgr.connect_sink("contacts", std::make_shared<dynamic::NdjsonFileSink>("contacts.jsonl"), sink_cfg);
        if (!json_out && cli.want_text) mgr.connect_sink("contacts", std::make_shared<dynamic::NdjsonFileSink>("contacts.txt"), sink_cfg);
        if (cli.want_log)  mgr.connect_sink("contacts", std::make_shared<dynamic::LoggingSink>(), sink_cfg);

        mgr.compile();
        mgr.run(static_cast<std::size_t>(cli.threads));
      }, src_variant);
    }

    Logger::get_logger()->info("Contact computation completed successfully!");
    return 0;

  } catch (const std::exception& e) {
    Logger::get_logger()->error("Error: {}", e.what());
    return 1;
  }
}

} // namespace lahuta::cli
