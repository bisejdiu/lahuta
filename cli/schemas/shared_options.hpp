/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto forward_concat = [](auto&& a, auto&& b, auto&& c) {
 *     return std::string(std::forward<decltype(a)>(a)) +
 *            std::forward<decltype(b)>(b) +
 *            std::forward<decltype(c)>(c);
 *   };
 *   return forward_concat("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_CLI_SHARED_OPTIONS_HPP
#define LAHUTA_CLI_SHARED_OPTIONS_HPP

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

#include "logging/logging.hpp"
#include "parsing/parsed_args.hpp"
#include "schemas/option_schema.hpp"

namespace lahuta::cli {

struct PipelineReporter;

namespace shared_opts {
constexpr unsigned BaseIndex = 100;
enum : unsigned {
  GlobalHelp = BaseIndex,
  GlobalVersion,
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
  bool version_requested             = false;
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

  RuntimeOptionSpec();
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
[[nodiscard]] std::filesystem::path validate_output_dir(const std::string &output_arg);
[[nodiscard]] std::string ensure_jsonl_extension(const std::string &path);
[[nodiscard]] std::string require_arg(const ParsedArgs &args,
                                      int option,
                                      std::string_view label,
                                      std::string_view missing_message = {},
                                      std::string_view empty_message = {});
[[nodiscard]] int default_thread_count();

} // namespace lahuta::cli

#endif // LAHUTA_CLI_SHARED_OPTIONS_HPP
