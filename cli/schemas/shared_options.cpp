#include <exception>
#include <sstream>
#include <string_view>

#include "parsing/arg_validation.hpp"
#include "parsing/extension_utils.hpp"
#include "parsing/usage_error.hpp"
#include "runner/reporting.hpp"
#include "schemas/shared_options.hpp"

namespace lahuta::cli {
namespace {

SourceOptionSpec default_source_spec() {
  SourceOptionSpec spec;
  spec.default_extensions = {".cif", ".cif.gz", ".pdb", ".pdb.gz"};
  return spec;
}

RuntimeOptionSpec default_runtime_spec() { return RuntimeOptionSpec{}; }

std::string join_option_list(const std::vector<std::string> &options) {
  std::ostringstream out;
  for (std::size_t i = 0; i < options.size(); ++i) {
    if (i > 0) {
      out << ", ";
    }
    out << options[i];
  }
  return out.str();
}

std::string build_extensions_help(const SourceOptionSpec &spec) {
  std::string help = "  --extension, -e <ext>        \tFile extension(s) for directory mode. Repeat or "
                     "comma-separate values";
  if (!spec.default_extensions.empty()) {
    help.append(" (default: ");
    help.append(describe_extensions(spec.default_extensions));
    help.append(").");
  } else {
    help.push_back('.');
  }
  return help;
}

std::string build_threads_help(int default_threads) {
  return "  --threads, -t <num>          \tNumber of threads to use (default: " +
         std::to_string(default_threads) + ").";
}

std::string build_batch_help(std::size_t default_batch_size) {
  return "  --batch-size, -b <size>      \tBatch size for processing (default: " +
         std::to_string(default_batch_size) + ").";
}

std::string build_writer_threads_help(std::size_t default_writer_threads) {
  return "  --writer-threads <num>       \tNumber of writer threads per sink (default: " +
         std::to_string(default_writer_threads) + ").";
}

std::size_t parse_size_t(std::string_view value, std::string_view label, bool allow_zero) {
  if (value.empty()) {
    throw CliUsageError(std::string(label) + " requires a value.");
  }
  std::size_t parsed = 0;
  try {
    parsed = std::stoull(std::string(value));
  } catch (const std::exception &) {
    throw CliUsageError("Invalid " + std::string(label) + " value '" + std::string(value) + "'");
  }
  if (!allow_zero && parsed == 0) {
    throw CliUsageError(std::string(label) + " must be positive");
  }
  return parsed;
}

int parse_int(std::string_view value, std::string_view label) {
  if (value.empty()) {
    throw CliUsageError(std::string(label) + " requires a value.");
  }
  int parsed = 0;
  try {
    parsed = std::stoi(std::string(value));
  } catch (const std::exception &) {
    throw CliUsageError("Invalid " + std::string(label) + " value '" + std::string(value) + "'");
  }
  if (parsed <= 0) {
    throw CliUsageError(std::string(label) + " must be positive");
  }
  return parsed;
}

} // namespace

void add_global_options(OptionSchema &schema) {
  schema.add({shared_opts::GlobalHelp,
              "h",
              "help",
              option::Arg::None,
              "  --help,  -h                  \tPrint this help message and exit."});
  schema.add({shared_opts::GlobalVersion,
              "",
              "version",
              option::Arg::None,
              "  --version                    \tPrint version information and exit."});
  schema.add({shared_opts::GlobalVerbose,
              "v",
              "verbose",
              validate::Verbosity,
              "  --verbose, -v <level>        \tSet verbosity level (0, 1, or 2)."});
  schema.add({shared_opts::GlobalProgressMs,
              "",
              "progress",
              validate::Required,
              "  --progress <ms>              \tProgress update in ms (0 disables, default: 50)."});
  schema.add({shared_opts::GlobalProgressNoColor,
              "",
              "progress-no-color",
              option::Arg::None,
              "  --progress-no-color          \tDisable progress output colors."});
}

void add_source_options(OptionSchema &schema) { add_source_options(schema, default_source_spec()); }

void add_source_options(OptionSchema &schema, const SourceOptionSpec &spec) {
  if (spec.allow_database) {
    schema.add({shared_opts::SourceDatabase,
                "",
                "database",
                validate::Required,
                "  --database <path>            \tProcess structures from database."});
  }
  if (spec.allow_directory) {
    schema.add({shared_opts::SourceDirectory,
                "d",
                "directory",
                validate::Required,
                "  --directory, -d <path>       \tProcess all files in directory."});
  }
  if (spec.allow_files) {
    schema.add({shared_opts::SourceVector,
                "f",
                "files",
                validate::Required,
                "  --files, -f <file1,file2>    \tProcess specific files (comma-separated or repeat -f)."});
  }
  if (spec.allow_file_list) {
    schema.add({shared_opts::SourceFileList,
                "l",
                "file-list",
                validate::Required,
                "  --file-list, -l <path>       \tProcess files listed in text file (one per line)."});
  }
  if (spec.allow_directory) {
    schema.add(
        {shared_opts::SourceExtension, "e", "extension", validate::Required, build_extensions_help(spec)});
    schema.add({shared_opts::SourceRecursive,
                "r",
                "recursive",
                option::Arg::None,
                "  --recursive, -r              \tRecursively search subdirectories."});
  }
  if (spec.include_is_af2_model) {
    schema.add({shared_opts::SourceIsAf2Model,
                "",
                "is_af2_model",
                option::Arg::None,
                "  --is_af2_model               \tInputs are AlphaFold2 models (or AF2-like)."});
  }
}

void add_runtime_options(OptionSchema &schema) { add_runtime_options(schema, default_runtime_spec()); }

void add_runtime_options(OptionSchema &schema, const RuntimeOptionSpec &spec) {
  schema.add({shared_opts::RuntimeThreads,
              "t",
              "threads",
              validate::Required,
              build_threads_help(spec.default_threads)});
  schema.add({shared_opts::RuntimeBatchSize,
              "b",
              "batch-size",
              validate::Required,
              build_batch_help(spec.default_batch_size)});
  if (spec.include_writer_threads) {
    schema.add({shared_opts::RuntimeWriterThreads,
                "",
                "writer-threads",
                validate::Required,
                build_writer_threads_help(spec.default_writer_threads)});
  }
}

void add_report_options(OptionSchema &schema) {
  schema.add({shared_opts::ReportReporter,
              "",
              "reporter",
              validate::Required,
              "  --reporter <name>            \tSelect pipeline reporter. See help for names."});
  schema.add({shared_opts::ReportSaveRunReport,
              "",
              "save-run-report",
              option::Arg::None,
              "  --save-run-report            \tSave run statistics to a JSON file."});
}

GlobalConfig parse_global_config(const ParsedArgs &args) {
  GlobalConfig config;
  config.help_requested    = args.get_flag(shared_opts::GlobalHelp);
  config.version_requested = args.get_flag(shared_opts::GlobalVersion);

  if (args.has(shared_opts::GlobalVerbose)) {
    const std::string level = args.get_string(shared_opts::GlobalVerbose);
    if (level == "0") {
      config.log_level = lahuta::Logger::LogLevel::Error;
    } else if (level == "1") {
      config.log_level = lahuta::Logger::LogLevel::Info;
    } else if (level == "2") {
      config.log_level = lahuta::Logger::LogLevel::Debug;
    } else if (level.empty()) {
      throw CliUsageError("Option '-v/--verbose' expects <level> (0,1,2)");
    } else {
      throw CliUsageError("Invalid verbosity level '" + level + "'. Must be 0, 1 or 2");
    }
  }

  if (args.has(shared_opts::GlobalProgressMs)) {
    config.progress_ms = parse_size_t(args.get_string(shared_opts::GlobalProgressMs), "--progress", true);
  }

  if (args.has(shared_opts::GlobalProgressNoColor)) {
    config.progress_color = false;
  }

  return config;
}

SourceConfig parse_source_config(const ParsedArgs &args) {
  return parse_source_config(args, default_source_spec());
}

SourceConfig parse_source_config(const ParsedArgs &args, const SourceOptionSpec &spec) {
  SourceConfig config;
  config.extensions = spec.default_extensions;

  int source_count = 0;
  if (spec.allow_database && args.has(shared_opts::SourceDatabase)) {
    config.mode          = SourceConfig::Mode::Database;
    config.database_path = args.get_string(shared_opts::SourceDatabase);
    source_count++;
  }
  if (spec.allow_directory && args.has(shared_opts::SourceDirectory)) {
    config.mode           = SourceConfig::Mode::Directory;
    config.directory_path = args.get_string(shared_opts::SourceDirectory);
    source_count++;
  }
  if (spec.allow_files && args.has(shared_opts::SourceVector)) {
    config.mode           = SourceConfig::Mode::Vector;
    const auto raw_values = args.get_all_strings(shared_opts::SourceVector);
    for (const auto &raw : raw_values) {
      parse_file_argument(raw, config.file_vector);
    }
    source_count++;
  }
  if (spec.allow_file_list && args.has(shared_opts::SourceFileList)) {
    config.mode           = SourceConfig::Mode::FileList;
    config.file_list_path = args.get_string(shared_opts::SourceFileList);
    source_count++;
  }

  if (source_count == 0) {
    std::vector<std::string> allowed;
    if (spec.allow_database) {
      allowed.emplace_back("--database");
    }
    if (spec.allow_directory) {
      allowed.emplace_back("--directory");
    }
    if (spec.allow_files) {
      allowed.emplace_back("--files");
    }
    if (spec.allow_file_list) {
      allowed.emplace_back("--file-list");
    }
    throw CliUsageError("Must specify exactly one source option: " + join_option_list(allowed));
  }
  if (source_count > 1) {
    throw CliUsageError("Cannot specify multiple source options");
  }

  if (spec.allow_directory && args.has(shared_opts::SourceExtension)) {
    config.extensions.clear();
    const auto raw_values = args.get_all_strings(shared_opts::SourceExtension);
    for (const auto &raw : raw_values) {
      parse_extension_argument(raw, config.extensions);
    }
    if (config.extensions.empty()) {
      config.extensions.emplace_back();
    }
  }

  if (spec.allow_directory) {
    config.recursive = args.get_flag(shared_opts::SourceRecursive);
  }

  if (spec.include_is_af2_model) {
    config.is_af2_model = args.get_flag(shared_opts::SourceIsAf2Model);
  }

  return config;
}

RuntimeConfig parse_runtime_config(const ParsedArgs &args) {
  return parse_runtime_config(args, default_runtime_spec());
}

RuntimeConfig parse_runtime_config(const ParsedArgs &args, const RuntimeOptionSpec &spec) {
  RuntimeConfig config;
  config.threads        = spec.default_threads;
  config.batch_size     = spec.default_batch_size;
  config.writer_threads = spec.default_writer_threads;

  if (args.has(shared_opts::RuntimeThreads)) {
    config.threads = parse_int(args.get_string(shared_opts::RuntimeThreads), "--threads");
  }

  if (args.has(shared_opts::RuntimeBatchSize)) {
    config.batch_size = parse_size_t(args.get_string(shared_opts::RuntimeBatchSize), "--batch-size", false);
  }

  if (spec.include_writer_threads && args.has(shared_opts::RuntimeWriterThreads)) {
    config.writer_threads = parse_size_t(args.get_string(shared_opts::RuntimeWriterThreads),
                                         "--writer-threads",
                                         false);
  }

  return config;
}

ReportConfig parse_report_config(const ParsedArgs &args) {
  ReportConfig config;
  config.reporter = &default_pipeline_reporter();

  if (args.has(shared_opts::ReportReporter)) {
    const std::string name = args.get_string(shared_opts::ReportReporter);
    if (name.empty()) {
      throw CliUsageError("--reporter requires a value.");
    }
    if (const auto *reporter = find_pipeline_reporter(name)) {
      config.reporter = reporter;
    } else {
      throw CliUsageError("Unknown reporter '" + name + "'");
    }
  }

  config.save_run_report = args.get_flag(shared_opts::ReportSaveRunReport);
  return config;
}

} // namespace lahuta::cli
