#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include "analysis/extract/extract_tasks.hpp"
#include "cli/arg_validation.hpp"
#include "cli/extension_utils.hpp"
#include "cli/run_report.hpp"
#include "cli/time_utils.hpp"
#include "commands/extract.hpp"
#include "commands/reporting.hpp"
#include "logging/logging.hpp"
#include "pipeline/dynamic/backpressure.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/dynamic/sources.hpp"
#include "runtime.hpp"
#include "sinks/logging.hpp"
#include "sinks/ndjson.hpp"

// clang-format off
namespace lahuta::cli {

namespace dyn = pipeline::dynamic;

struct ExtractOptions {
  enum class SourceMode { Directory, Vector, FileList, Database };

  SourceMode source_mode = SourceMode::Directory;
  std::string directory_path;
  std::vector<std::string> extensions{".cif", ".cif.gz", ".pdb", ".pdb.gz"};
  bool recursive = false;
  std::vector<std::string> file_vector;
  std::string file_list_path;
  std::string database_path;
  std::vector<std::string> fields;

  bool is_af2_model = false;
  bool output_stdout = false;
  std::string output_path;
  bool output_override = false;
  bool save_run_report = false;
  const PipelineReporter* reporter = nullptr;
  int threads = 8;
  size_t batch_size = 512;
  size_t writer_threads = 1;
};

namespace {

constexpr std::string_view FIELD_SEQUENCE = "sequence";
constexpr std::string_view FIELD_PLDDT    = "plddt";
constexpr std::string_view FIELD_DSSP     = "dssp";
constexpr std::string_view FIELD_ORGANISM = "organism";

bool is_valid_field(std::string_view field) {
  return field == FIELD_SEQUENCE ||
         field == FIELD_PLDDT ||
         field == FIELD_DSSP ||
         field == FIELD_ORGANISM;
}

void parse_fields_argument(std::string_view raw, std::vector<std::string>& out) {
  detail::split_argument_list(raw, false, out);
}

void initialize_runtime(int num_threads) {
  LahutaRuntime::ensure_initialized(static_cast<std::size_t>(num_threads));
}

std::shared_ptr<dyn::ITask> build_extract_task(std::string_view field, std::string output_channel) {
  using namespace analysis::extract;
  if (field == FIELD_SEQUENCE) {
    return std::make_shared<SequenceExtractTask>(std::move(output_channel));
  }
  if (field == FIELD_PLDDT) {
    return std::make_shared<PlddtExtractTask>(std::move(output_channel));
  }
  if (field == FIELD_DSSP) {
    return std::make_shared<DsspExtractTask>(std::move(output_channel));
  }
  if (field == FIELD_ORGANISM) {
    return std::make_shared<OrganismExtractTask>(std::move(output_channel));
  }
  throw std::logic_error("Invalid extract field");
}

} // namespace

namespace extract_opts {
const option::Descriptor usage[] = {
  {ExtractOptionIndex::Unknown, 0, "", "", validate::Unknown,
   "Usage: lahuta extract --fields <fields> [options]\n\n"
   "Extract data from AlphaFold2 model files or databases.\n"
   "Note: file-based inputs must be AF2 model files (mmCIF). Generic structures are not supported.\n\n"
   "Available reporters (--reporter <name>):\n"
   "  summary     - Balanced totals and item counts (default; negligible overhead).\n"
   "  terse       - Single-line throughput summary (fastest logging footprint).\n"
   "  diagnostics - Summary plus concurrency and stage breakdown details.\n\n"
   "Fields (repeat --fields or comma-separate):\n"
   "  sequence   Extract amino acid sequences.\n"
   "  plddt      Extract per-residue pLDDT confidence scores.\n"
   "  dssp       Extract secondary structure assignments (DSSP).\n"
   "  organism   Extract organism metadata.\n\n"
   "Input Options (choose one):"},
  {ExtractOptionIndex::Help, 0, "h", "help", option::Arg::None,
   "  --help, -h                   \tPrint this help message and exit."},
  {ExtractOptionIndex::Fields, 0, "", "fields", validate::Required,
   "  --fields <field1,field2>     \tFields to extract (repeat or comma-separate)."},
  {ExtractOptionIndex::SourceDatabase, 0, "", "database", validate::Required,
   "  --database <path>            \tProcess structures from database."},
  {ExtractOptionIndex::SourceDirectory, 0, "d", "directory", validate::Required,
   "  --directory, -d <path>       \tProcess all files in directory."},
  {ExtractOptionIndex::SourceVector, 0, "f", "files", validate::Required,
   "  --files, -f <file1,file2>    \tProcess specific files (comma-separated or repeat -f)."},
  {ExtractOptionIndex::SourceFileList, 0, "l", "file-list", validate::Required,
   "  --file-list, -l <path>       \tProcess files listed in text file (one per line)."},
  {0, 0, "", "", option::Arg::None,
   "\nDirectory Options:"},
  {ExtractOptionIndex::Extension, 0, "e", "extension", validate::Required,
   "  --extension, -e <ext>        \tFile extension(s) for directory mode. Repeat or comma-separate values (default: .cif, .cif.gz, .pdb, .pdb.gz)."},
  {ExtractOptionIndex::Recursive, 0, "r", "recursive", option::Arg::None,
   "  --recursive, -r              \tRecursively search subdirectories."},
  {0, 0, "", "", option::Arg::None,
   "\nModel Options:"},
  {ExtractOptionIndex::IsAf2Model, 0, "", "is_af2_model", option::Arg::None,
   "  --is_af2_model               \tRequired for file-based sources; inputs are AlphaFold2 models (AF2-like mmCIF)."},
  {0, 0, "", "", option::Arg::None,
   "\nOutput Options:"},
  {ExtractOptionIndex::Output, 0, "o", "output", validate::Required,
   "  --output, -o <path>          \tWrite NDJSON to file (default: <field>_data.jsonl). Use '-' for stdout."},
  {ExtractOptionIndex::Reporter, 0, "", "reporter", validate::Required,
   "  --reporter <name>            \tSelect pipeline reporter. See help for names."},
  {ExtractOptionIndex::SaveRunReport, 0, "", "save-run-report", option::Arg::None,
   "  --save-run-report            \tSave run statistics to a JSON file."},
  {0, 0, "", "", option::Arg::None,
   "\nRuntime Options:"},
  {ExtractOptionIndex::Threads, 0, "t", "threads", validate::Required,
   "  --threads, -t <num>          \tNumber of threads to use (default: 8)."},
  {ExtractOptionIndex::BatchSize, 0, "b", "batch-size", validate::Required,
   "  --batch-size, -b <size>      \tBatch size for processing (default: 512)."},
  {ExtractOptionIndex::WriterThreads, 0, "", "writer-threads", validate::Required,
   "  --writer-threads <num>       \tNumber of writer threads per sink (default: 1)."},
  {0, 0, 0, 0, 0, 0}
};
} // namespace extract_opts

[[nodiscard]] std::unique_ptr<CliCommand> ExtractCommand::create() {
  return std::unique_ptr<CliCommand>(new ExtractCommand());
}

int ExtractCommand::run(int argc, char* argv[]) {
  option::Stats stats(true, extract_opts::usage, argc, const_cast<const char**>(argv));
  std::vector<option::Option> options(stats.options_max);
  std::vector<option::Option> buffer (stats.buffer_max);
  option::Parser parse(true, extract_opts::usage, argc, const_cast<const char**>(argv), options.data(), buffer.data());

  if (parse.error()) return 1;

  if (options[extract_opts::ExtractOptionIndex::Help]) {
    option::printUsage(std::cout, extract_opts::usage);
    return 0;
  }

  try {
    if (parse.nonOptionsCount() > 0) {
      Logger::get_logger()->error("Unexpected positional argument '{}'. Use --fields instead.", parse.nonOption(0));
      option::printUsage(std::cerr, extract_opts::usage);
      return 1;
    }

    std::vector<std::string> field_tokens;
    for (const option::Option* opt = &options[extract_opts::ExtractOptionIndex::Fields];
         opt != nullptr;
         opt = opt->next()) {
      parse_fields_argument(opt->arg ? opt->arg : "", field_tokens);
    }

    if (field_tokens.empty()) {
      Logger::get_logger()->error("Missing required --fields option. Must include one or more of: sequence, plddt, dssp, organism");
      option::printUsage(std::cerr, extract_opts::usage);
      return 1;
    }

    std::vector<std::string> fields;
    std::unordered_set<std::string> seen;
    for (const auto& token : field_tokens) {
      if (!is_valid_field(token)) {
        Logger::get_logger()->error("Invalid field '{}'. Must be one of: sequence, plddt, dssp, organism", token);
        option::printUsage(std::cerr, extract_opts::usage);
        return 1;
      }
      if (seen.insert(token).second) {
        fields.push_back(token);
      }
    }
    const std::string run_timestamp = current_timestamp_string();

    ExtractOptions cli;
    const auto default_sink_cfg = dyn::get_default_backpressure_config();
    cli.writer_threads = default_sink_cfg.writer_threads;
    cli.reporter = &default_pipeline_reporter();
    cli.output_path.clear();
    cli.output_override = false;
    cli.fields = std::move(fields);

    // Parse source options
    int source_count = 0;
    if (options[extract_opts::ExtractOptionIndex::SourceDatabase]) {
      cli.source_mode = ExtractOptions::SourceMode::Database;
      cli.database_path = options[extract_opts::ExtractOptionIndex::SourceDatabase].arg;
      source_count++;
    }
    if (options[extract_opts::ExtractOptionIndex::SourceDirectory]) {
      cli.source_mode = ExtractOptions::SourceMode::Directory;
      cli.directory_path = options[extract_opts::ExtractOptionIndex::SourceDirectory].arg;
      source_count++;
    }
    if (options[extract_opts::ExtractOptionIndex::SourceVector]) {
      cli.source_mode = ExtractOptions::SourceMode::Vector;
      cli.file_vector.clear();
      for (const option::Option* opt = &options[extract_opts::ExtractOptionIndex::SourceVector];
           opt != nullptr;
           opt = opt->next()) {
        if (opt->arg) parse_file_argument(opt->arg, cli.file_vector);
      }
      source_count++;
    }
    if (options[extract_opts::ExtractOptionIndex::SourceFileList]) {
      cli.source_mode = ExtractOptions::SourceMode::FileList;
      cli.file_list_path = options[extract_opts::ExtractOptionIndex::SourceFileList].arg;
      source_count++;
    }

    if (source_count == 0) {
      Logger::get_logger()->error("Must specify exactly one source option: --database, --directory, --files, or --file-list");
      return 1;
    }
    if (source_count > 1) {
      Logger::get_logger()->error("Cannot specify multiple source options");
      return 1;
    }

    if (options[extract_opts::ExtractOptionIndex::Extension]) {
      cli.extensions.clear();
      for (const option::Option* opt = &options[extract_opts::ExtractOptionIndex::Extension]; opt != nullptr; opt = opt->next()) {
        parse_extension_argument(opt->arg ? opt->arg : "", cli.extensions);
      }
    }
    if (cli.extensions.empty()) cli.extensions.emplace_back();

    cli.recursive    = options[extract_opts::ExtractOptionIndex::Recursive]  ? true : false;
    cli.is_af2_model = options[extract_opts::ExtractOptionIndex::IsAf2Model] ? true : false;

    if (cli.source_mode == ExtractOptions::SourceMode::Database) {
      cli.is_af2_model = true;
    }

    if (cli.source_mode != ExtractOptions::SourceMode::Database && !cli.is_af2_model) {
      Logger::get_logger()->error(
          "extract expects AlphaFold2 model inputs. For file-based sources, pass --is_af2_model (or use --database).");
      option::printUsage(std::cerr, extract_opts::usage);
      return 1;
    }

    if (options[extract_opts::ExtractOptionIndex::Threads]) {
      cli.threads = std::stoi(options[extract_opts::ExtractOptionIndex::Threads].arg);
      if (cli.threads <= 0) {
        Logger::get_logger()->error("Threads must be positive");
        return 1;
      }
    }

    if (options[extract_opts::ExtractOptionIndex::BatchSize]) {
      cli.batch_size = std::stoull(options[extract_opts::ExtractOptionIndex::BatchSize].arg);
      if (cli.batch_size == 0) {
        Logger::get_logger()->error("Batch size must be positive");
        return 1;
      }
    }

    if (options[extract_opts::ExtractOptionIndex::WriterThreads]) {
      cli.writer_threads = std::stoull(options[extract_opts::ExtractOptionIndex::WriterThreads].arg);
      if (cli.writer_threads == 0) {
        Logger::get_logger()->error("Writer threads must be positive");
        return 1;
      }
    }

    if (options[extract_opts::ExtractOptionIndex::Output]) {
      std::string_view output_arg = options[extract_opts::ExtractOptionIndex::Output].arg
                                      ? options[extract_opts::ExtractOptionIndex::Output].arg
                                      : std::string_view{};
      if (output_arg.empty()) {
        Logger::get_logger()->error("--output requires a value.");
        return 1;
      }
      if (output_arg == "-") {
        cli.output_stdout = true;
        cli.output_path = "-";
      } else {
        cli.output_path = std::string(output_arg);
        cli.output_override = true;
      }
    }

    if (options[extract_opts::ExtractOptionIndex::Reporter]) {
      std::string_view name = options[extract_opts::ExtractOptionIndex::Reporter].arg
                                ? options[extract_opts::ExtractOptionIndex::Reporter].arg
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

    cli.save_run_report = options[extract_opts::ExtractOptionIndex::SaveRunReport] ? true : false;

    if (cli.fields.size() > 1 && cli.output_override) {
      Logger::get_logger()->error("--output can only be used with a single field. Omit --output for per-field outputs or use --output - for stdout.");
      return 1;
    }

    initialize_runtime(cli.threads);

    auto source = std::unique_ptr<sources::IDescriptor>{};
    switch (cli.source_mode) {
      case ExtractOptions::SourceMode::Database:
        source = dyn::sources_factory::from_lmdb(
            cli.database_path,
            std::string{},
            cli.batch_size,
            {static_cast<std::size_t>(cli.threads) + 1});
        break;
      case ExtractOptions::SourceMode::Directory:
        source = dyn::sources_factory::from_directory(
            cli.directory_path,
            cli.extensions,
            cli.recursive,
            cli.batch_size);
        break;
      case ExtractOptions::SourceMode::Vector:
        source = dyn::sources_factory::from_vector(cli.file_vector);
        break;
      case ExtractOptions::SourceMode::FileList:
        source = dyn::sources_factory::from_filelist(cli.file_list_path);
        break;
    }

    dyn::StageManager mgr(std::move(source));
    mgr.set_reporting_level(reporting_level_for_reporter(cli.reporter));

    auto sink_cfg = dyn::get_default_backpressure_config();
    sink_cfg.writer_threads = cli.writer_threads;

    const bool needs_parse_task = cli.source_mode != ExtractOptions::SourceMode::Database;
    if (needs_parse_task) {
      auto parse_task = std::make_shared<analysis::extract::ModelParseTask>();
      mgr.add_task("parse_model", {}, std::move(parse_task), /*thread_safe=*/true);
    }

    std::shared_ptr<dyn::LoggingSink> stdout_sink;
    if (cli.output_stdout) {
      stdout_sink = std::make_shared<dyn::LoggingSink>();
      Logger::get_logger()->info("Extract output sink -> stdout (-)");
    }

    for (const auto& field : cli.fields) {
      const std::string output_channel = field + "_data";
      auto task = build_extract_task(field, output_channel);
      std::vector<std::string> deps;
      if (needs_parse_task) deps.emplace_back("parse_model");
      mgr.add_task(output_channel, std::move(deps), std::move(task), /*thread_safe=*/true);

      if (cli.output_stdout) {
        mgr.connect_sink(output_channel, stdout_sink, sink_cfg);
      } else {
        const std::string output_path = cli.output_override ? cli.output_path : output_channel + ".jsonl";
        mgr.connect_sink(output_channel, std::make_shared<dyn::NdjsonFileSink>(output_path), sink_cfg);
        Logger::get_logger()->info("Extract output sink -> file: {}", output_path);
      }
    }

    mgr.compile();
    auto progress = attach_progress_observer(mgr, "extract");
    const auto report = mgr.run(static_cast<std::size_t>(cli.threads));
    if (progress) progress->finish();
    const auto* reporter = cli.reporter ? cli.reporter : &default_pipeline_reporter();
    reporter->emit("extract", report);
    if (cli.save_run_report) {
      const std::string report_path = make_report_path("extract", report.run_token, run_timestamp);
      if (!write_run_report_json(report_path, report)) {
        throw std::runtime_error("Failed to persist RunReport JSON");
      }
      Logger::get_logger()->info("Run report saved to {}", report_path);
    }
    return 0;
  } catch (const std::exception& e) {
    Logger::get_logger()->error("Error: {}", e.what());
    return 1;
  }
}

} // namespace lahuta::cli
