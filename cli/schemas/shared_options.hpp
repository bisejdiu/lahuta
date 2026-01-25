#ifndef LAHUTA_CLI_SHARED_OPTIONS_HPP
#define LAHUTA_CLI_SHARED_OPTIONS_HPP

#include <cstddef>
#include <string>
#include <vector>

#include "parsing/parsed_args.hpp"
#include "schemas/option_schema.hpp"
#include "logging/logging.hpp"

namespace lahuta::cli {

struct PipelineReporter;

namespace shared_opts {
constexpr unsigned BaseIndex = 100;
enum OptionIndex : unsigned {
  GlobalHelp = BaseIndex,
  GlobalVerbose,
  GlobalProgressMs,
  GlobalProgressNoColor,
  SourceDatabase,
  SourceDirectory,
  SourceVector,
  SourceFileList,
  SourceExtension,
  SourceRecursive,
  SourceIsAf2Model,
  RuntimeThreads,
  RuntimeBatchSize,
  RuntimeWriterThreads,
  ReportReporter,
  ReportSaveRunReport
};
} // namespace shared_opts

struct GlobalConfig {
  lahuta::Logger::LogLevel log_level = lahuta::Logger::LogLevel::Info;
  std::size_t progress_ms            = 50;
  bool progress_color                = true;
  bool help_requested                = false;
};

struct SourceOptionSpec {
  bool allow_database       = true;
  bool allow_directory      = true;
  bool allow_files          = true;
  bool allow_file_list      = true;
  bool include_is_af2_model = false;
  std::vector<std::string> default_extensions;
};

struct SourceConfig {
  enum class Mode { Directory, Vector, FileList, Database };

  Mode mode = Mode::Directory;
  std::string directory_path;
  std::vector<std::string> extensions;
  bool recursive = false;
  std::vector<std::string> file_vector;
  std::string file_list_path;
  std::string database_path;
  bool is_af2_model = false;
};

struct RuntimeOptionSpec {
  bool include_writer_threads        = false;
  int default_threads                = 8;
  std::size_t default_batch_size     = 512;
  std::size_t default_writer_threads = 1;
};

struct RuntimeConfig {
  int threads                = 8;
  std::size_t batch_size     = 512;
  std::size_t writer_threads = 1;
};

struct ReportConfig {
  const PipelineReporter *reporter = nullptr;
  bool save_run_report             = false;
};

void add_global_options(OptionSchema &schema);
void add_source_options(OptionSchema &schema);
void add_source_options(OptionSchema &schema, const SourceOptionSpec &spec);
void add_runtime_options(OptionSchema &schema);
void add_runtime_options(OptionSchema &schema, const RuntimeOptionSpec &spec);
void add_report_options(OptionSchema &schema);

[[nodiscard]] GlobalConfig parse_global_config(const ParsedArgs &args);
[[nodiscard]] SourceConfig parse_source_config(const ParsedArgs &args);
[[nodiscard]] SourceConfig parse_source_config(const ParsedArgs &args, const SourceOptionSpec &spec);
[[nodiscard]] RuntimeConfig parse_runtime_config(const ParsedArgs &args);
[[nodiscard]] RuntimeConfig parse_runtime_config(const ParsedArgs &args, const RuntimeOptionSpec &spec);
[[nodiscard]] ReportConfig parse_report_config(const ParsedArgs &args);

} // namespace lahuta::cli

#endif // LAHUTA_CLI_SHARED_OPTIONS_HPP
