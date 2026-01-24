#include <algorithm>
#include <filesystem>
#include <iostream>
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
#include "cli/run_report.hpp"
#include "cli/time_utils.hpp"
#include "commands/contacts.hpp"
#include "commands/reporting.hpp"
#include "db/db.hpp"
#include "logging/logging.hpp"
#include "pipeline/compute/parameters.hpp"
#include "pipeline/dynamic/backpressure.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/dynamic/sources.hpp"
#include "runtime.hpp"
#include "sinks/logging.hpp"
#include "sinks/ndjson.hpp"

// clang-format off
namespace lahuta::cli {

using namespace lahuta::pipeline;

using Source = std::variant<sources::Directory, std::vector<std::string>, sources::FileList>;

static void initialize_runtime(int num_threads) {
  lahuta::LahutaRuntime::ensure_initialized(static_cast<std::size_t>(num_threads));
}

struct ContactsOptions {
  enum class SourceMode { Directory, Vector, FileList, Database, MD };

  SourceMode source_mode = SourceMode::Directory;
  std::string directory_path;
  std::vector<std::string> extensions{".cif", ".cif.gz"};
  bool recursive = true;
  std::vector<std::string> file_vector;
  std::string file_list_path;
  std::string database_path;

  // MD fields
  std::string md_structure_path;                // PDB or GRO file
  std::vector<std::string> md_trajectory_paths; // XTC files

  analysis::contacts::ContactProvider provider = analysis::contacts::ContactProvider::MolStar;
  InteractionTypeSet interaction_types = InteractionTypeSet::all();

  bool is_af2_model = false;
  bool want_json = false;
  bool want_text = false;
  bool want_log  = false;
  bool save_run_report = false;
  int threads = 8;
  size_t batch_size = 200;
  size_t writer_threads = 1;
  const PipelineReporter* reporter = nullptr;
};

Source pick_source(const ContactsOptions& cli) {
  switch (cli.source_mode) {
    case ContactsOptions::SourceMode::Directory: return sources::Directory{cli.directory_path, cli.extensions, cli.recursive, cli.batch_size};
    case ContactsOptions::SourceMode::Vector:    return cli.file_vector;
    case ContactsOptions::SourceMode::FileList:  return sources::FileList {cli.file_list_path};
    case ContactsOptions::SourceMode::Database:  break; // Handled separately
    case ContactsOptions::SourceMode::MD:        break; // Handled separately
  }
  throw std::logic_error("Invalid source mode");
}

bool has_md_extension(const std::string& path, const std::vector<std::string>& extensions) { // case-insensitive
  std::filesystem::path p(path);
  std::string ext = p.extension().string();
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
  for (const auto& valid_ext : extensions) {
    if (ext == valid_ext) return true;
  }
  return false;
}

// MD trajectory source descriptor
// produces a single IngestDescriptor with MDRef
class SingleTrajectoryDescriptor final : public lahuta::sources::IDescriptor {
public:
  SingleTrajectoryDescriptor(std::string structure, std::vector<std::string> xtcs)
      : structure_(std::move(structure)), xtcs_(std::move(xtcs)) {}

  std::optional<IngestDescriptor> next() override {
    if (done_) return std::nullopt;
    done_ = true;
    IngestDescriptor desc;
    desc.id = structure_;
    desc.origin = MDRef{structure_, xtcs_};
    return desc;
  }

  void reset() override { done_ = false; }

private:
  std::string structure_;
  std::vector<std::string> xtcs_;
  bool done_ = false;
};

namespace contacts_opts {
const option::Descriptor usage[] = {
  {ContactsOptionIndex::Unknown, 0, "", "", validate::Unknown,
   "Usage: lahuta contacts [options]\n"
   "Author: Besian I. Sejdiu (@bisejdiu)\n\n"
   "Compute inter-atomic contacts.\n\n"

    "<<< Detailed help will be added soon. >>>\n\n"

    "Available reporters (--reporter <name>):\n"
    "  summary     - Concise summary of computation results (negligible overhead).\n"
    "  terse       - Fastest logging footprint.\n"
    "  diagnostics - Detailed diagnostics for in-depth analysis (small overhead).\n"
    "\n"

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
  {ContactsOptionIndex::SourceMD, 0, "", "md", validate::Required,
   "  --md <path>                  \tMD input file. Specify once for structure (PDB/GRO) and once or more for trajectories (XTC)."},
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
  {ContactsOptionIndex::IsAf2Model, 0, "", "is_af2_model", option::Arg::None,
   "  --is_af2_model               \tInputs are AlphaFold2 models (or AF2-like)."},
  {0, 0, "", "", option::Arg::None,
   "\nOutput Options:"},
  {ContactsOptionIndex::OutputJson, 0, "", "json", option::Arg::None,
   "  --json                       \tOutput results in JSON format."},
  {ContactsOptionIndex::OutputText, 0, "", "text", option::Arg::None,
   "  --text                       \tOutput results in text format."},
  {ContactsOptionIndex::OutputLog, 0, "", "log", option::Arg::None,
   "  --log                        \tOutput results to standard output (logging)."},
  {ContactsOptionIndex::Reporter, 0, "", "reporter", validate::Required,
   "  --reporter <name>            \tSelect pipeline reporter. See help for names."},
  {ContactsOptionIndex::SaveRunReport, 0, "", "save-run-report", option::Arg::None,
   "  --save-run-report            \tSave run statistics to a JSON file."},
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
    const auto default_sink_cfg = dynamic::get_default_backpressure_config();
    cli.writer_threads = default_sink_cfg.writer_threads;
    cli.reporter = &default_pipeline_reporter();
    const std::string run_timestamp = current_timestamp_string();
    const std::string contacts_json_path = "contacts_" + run_timestamp + ".jsonl";
    const std::string contacts_text_path = "contacts_" + run_timestamp + ".txt";

    auto emit_and_save_report = [&](const pipeline::dynamic::StageManager::RunReport& report) {
      const auto* reporter = cli.reporter ? cli.reporter : &default_pipeline_reporter();
      reporter->emit("contacts", report);
      if (cli.save_run_report) {
        const std::string report_path = make_report_path("contacts", report.run_token, run_timestamp);
        if (!write_run_report_json(report_path, report)) {
          throw std::runtime_error("Failed to persist RunReport JSON");
        }
        Logger::get_logger()->info("Run report saved to {}", report_path);
      }
    };

    if (options[contacts_opts::ContactsOptionIndex::Reporter]) {
      std::string_view name = options[contacts_opts::ContactsOptionIndex::Reporter].arg
                                ? options[contacts_opts::ContactsOptionIndex::Reporter].arg
                                : std::string_view{};
      if (name.empty()) {
        Logger::get_logger()->error("--reporter requires a value.");
        return 1;
      }
      if (const auto* rep = find_pipeline_reporter(name)) {
        cli.reporter = rep;
      } else {
        Logger::get_logger()->error("Unknown reporter '{}' ", name);
        return 1;
      }
    }

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
    if (options[contacts_opts::ContactsOptionIndex::SourceMD]) {
      cli.source_mode = ContactsOptions::SourceMode::MD;
      // Collect all --md arguments
      std::vector<std::string> md_files;
      for (const option::Option* opt = &options[contacts_opts::ContactsOptionIndex::SourceMD];
           opt != nullptr;
           opt = opt->next()) {
        if (opt->arg) md_files.push_back(opt->arg);
      }

      // Validate and separate structure vs. trajectory files
      std::vector<std::string> structures;
      std::vector<std::string> trajectories;
      for (const auto& file : md_files) {
        if (!std::filesystem::exists(file)) {
          Logger::get_logger()->error("MD input file not found: {}", file);
          return 1;
        }
        if (has_md_extension(file, {".pdb", ".gro"})) {
          structures.push_back(file);
        } else if (has_md_extension(file, {".xtc"})) {
          trajectories.push_back(file);
        } else {
          Logger::get_logger()->error("Invalid MD file extension: {}. Expected .pdb, .gro, or .xtc", file);
          return 1;
        }
      }

      // Validate: exactly one structure file required
      if (structures.empty()) {
        Logger::get_logger()->error("MD mode requires exactly one structure file (.pdb or .gro)");
        return 1;
      }
      if (structures.size() > 1) {
        Logger::get_logger()->error("MD mode requires exactly one structure file, but {} were provided", structures.size());
        return 1;
      }

      if (trajectories.empty()) {
        Logger::get_logger()->error("MD mode requires at least one trajectory file (.xtc)");
        return 1;
      }

      cli.md_structure_path = structures[0];
      cli.md_trajectory_paths = std::move(trajectories);
      source_count++;
    }

    if (source_count == 0) {
      Logger::get_logger()->error("Must specify exactly one source option: --directory, --files, --file-list, --database, or --md");
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

    cli.is_af2_model = options[contacts_opts::ContactsOptionIndex::IsAf2Model] ? true : false;

    // Parse output options
    cli.want_json        = options[contacts_opts::ContactsOptionIndex::OutputJson]     ? true : false;
    cli.want_text        = options[contacts_opts::ContactsOptionIndex::OutputText]     ? true : false;
    cli.want_log         = options[contacts_opts::ContactsOptionIndex::OutputLog]      ? true : false;
    cli.save_run_report  = options[contacts_opts::ContactsOptionIndex::SaveRunReport]  ? true : false;

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
    bool is_md = (cli.source_mode == ContactsOptions::SourceMode::MD);
    const bool use_model_pipeline = cli.is_af2_model || is_db;

    if (is_md) {
      Logger::get_logger()->info("MD Structure file: {}", cli.md_structure_path);
      Logger::get_logger()->info("MD Trajectory files ({}):", cli.md_trajectory_paths.size());
      for (const auto& xtc : cli.md_trajectory_paths) {
        Logger::get_logger()->info("  - {}", xtc);
      }

      auto source = std::make_unique<SingleTrajectoryDescriptor>(cli.md_structure_path, cli.md_trajectory_paths);
      dynamic::StageManager mgr(std::move(source));
      mgr.set_auto_builtins(true);
      mgr.set_reporting_level(reporting_level_for_reporter(cli.reporter));

      // Configure built-ins
      mgr.get_topology_params().atom_typing_method = analysis::contacts::typing_for_provider(cli.provider);

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
      auto sink_cfg = dynamic::get_default_backpressure_config();
      sink_cfg.writer_threads = cli.writer_threads;
      if  (json_out && cli.want_json) {
        mgr.connect_sink("contacts", std::make_shared<dynamic::NdjsonFileSink>(contacts_json_path), sink_cfg);
        Logger::get_logger()->info("Contacts JSON sink -> {}", contacts_json_path);
      }
      if (!json_out && cli.want_text) {
        mgr.connect_sink("contacts", std::make_shared<dynamic::NdjsonFileSink>(contacts_text_path), sink_cfg);
        Logger::get_logger()->info("Contacts text sink -> {}", contacts_text_path);
      }
      if (cli.want_log) mgr.connect_sink("contacts", std::make_shared<dynamic::LoggingSink>(), sink_cfg);

      mgr.compile();
      const auto report = mgr.run(static_cast<std::size_t>(cli.threads));
      emit_and_save_report(report);
    } else if (is_db) {
      auto src = dynamic::sources_factory::from_lmdb(
          cli.database_path,
          std::string{},
          cli.batch_size,
          {static_cast<std::size_t>(cli.threads) + 1}); // +1 for the main thread reading keys.
      dynamic::StageManager mgr(std::move(src));
      // Enable conditional built-ins injection for CLI ergonomics
      mgr.set_auto_builtins(true);
      mgr.set_reporting_level(reporting_level_for_reporter(cli.reporter));

      // Configure built-ins
      mgr.get_system_params().is_model = use_model_pipeline;
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
      auto sink_cfg = dynamic::get_default_backpressure_config();
      sink_cfg.writer_threads = cli.writer_threads;
      if  (json_out && cli.want_json) {
        mgr.connect_sink("contacts", std::make_shared<dynamic::NdjsonFileSink>(contacts_json_path), sink_cfg);
        Logger::get_logger()->info("Contacts JSON sink -> {}", contacts_json_path);
      }
      if (!json_out && cli.want_text) {
        mgr.connect_sink("contacts", std::make_shared<dynamic::NdjsonFileSink>(contacts_text_path), sink_cfg);
        Logger::get_logger()->info("Contacts text sink -> {}", contacts_text_path);
      }
      if (cli.want_log)  mgr.connect_sink("contacts", std::make_shared<dynamic::LoggingSink>(), sink_cfg);

      mgr.compile();
      auto progress = attach_progress_observer(mgr);
      const auto report = mgr.run(static_cast<std::size_t>(cli.threads));
      if (progress) progress->finish();
      emit_and_save_report(report);
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
        mgr.set_reporting_level(reporting_level_for_reporter(cli.reporter));

        // Configure built-ins
        mgr.get_system_params().is_model = use_model_pipeline;
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
        auto sink_cfg = dynamic::get_default_backpressure_config();
        sink_cfg.writer_threads = cli.writer_threads;
        if  (json_out && cli.want_json) {
          mgr.connect_sink("contacts", std::make_shared<dynamic::NdjsonFileSink>(contacts_json_path), sink_cfg);
          Logger::get_logger()->info("Contacts JSON sink -> {}", contacts_json_path);
        }
        if (!json_out && cli.want_text) {
          mgr.connect_sink("contacts", std::make_shared<dynamic::NdjsonFileSink>(contacts_text_path), sink_cfg);
          Logger::get_logger()->info("Contacts text sink -> {}", contacts_text_path);
        }
        if (cli.want_log)  mgr.connect_sink("contacts", std::make_shared<dynamic::LoggingSink>(), sink_cfg);

        mgr.compile();
        auto progress = attach_progress_observer(mgr);
        const auto report = mgr.run(static_cast<std::size_t>(cli.threads));
        if (progress) progress->finish();
        emit_and_save_report(report);
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
