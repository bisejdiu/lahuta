/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto make = []() noexcept(noexcept(std::string{})) {
 *     return std::string{"besian"} + "sejdiu" + "@gmail.com";
 *   };
 *   static_assert(noexcept(make()) == noexcept(std::string{}));
 *   return make();
 * }();
 *
 */

#include <memory>
#include <string>
#include <string_view>
#include <utility>

#include "analysis/system/model_pack_task.hpp"
#include "db/db.hpp"
#include "logging/logging.hpp"
#include "parsing/arg_validation.hpp"
#include "parsing/extension_utils.hpp"
#include "parsing/usage_error.hpp"
#include "pipeline/ingest/factory.hpp"
#include "runner/time_utils.hpp"
#include "schemas/shared_options.hpp"
#include "sinks/lmdb.hpp"
#include "specs/command_spec.hpp"

namespace lahuta::cli {
namespace A = lahuta::analysis;
namespace P = lahuta::pipeline;
namespace {

constexpr std::string_view Summary = "Create a database from alphafold2 models.";

namespace createdb_opts {
constexpr unsigned BaseIndex = 200;
enum : unsigned { OutputPath = BaseIndex, MaxSize };
} // namespace createdb_opts

constexpr std::size_t DefaultMaxSizeGb = 500;

struct CreateDbConfig {
  SourceConfig source;
  RuntimeConfig runtime;
  ReportConfig report;
  std::string database_path;
  std::size_t max_size_gb = DefaultMaxSizeGb;
};

std::size_t parse_size_option(std::string_view value, std::string_view label) {
  if (value.empty()) {
    throw CliUsageError(std::string(label) + " requires a value.");
  }
  std::size_t parsed = 0;
  try {
    parsed = std::stoull(std::string(value));
  } catch (const std::exception &) {
    throw CliUsageError("Invalid " + std::string(label) + " value '" + std::string(value) + "'");
  }
  if (parsed == 0) {
    throw CliUsageError(std::string(label) + " must be positive");
  }
  return parsed;
}

class CreateDbSpec final : public CommandSpec {
public:
  CreateDbSpec() {
    source_spec_.allow_database       = false;
    source_spec_.include_is_af2_model = false;
    source_spec_.default_extensions   = {".cif", ".cif.gz"};

    runtime_spec_.include_writer_threads = false;
    runtime_spec_.default_batch_size     = 1000;

    schema_.add({0,
                 "",
                 "",
                 validate::Unknown,
                 std::string("Usage: lahuta createdb [options]\n"
                             "Author: ")
                     .append(Author)
                     .append("\n\n")
                     .append(Summary)});

    schema_.add({0, "", "", option::Arg::None, "\nInput Options (choose one):"});
    add_source_options(schema_, source_spec_);

    schema_.add({0, "", "", option::Arg::None, "\nOutput Options:"});
    schema_.add({createdb_opts::OutputPath,
                 "o",
                 "output",
                 validate::Required,
                 "  --output, -o <path>          \tOutput database path "
                 "(default: createdb_<timestamp>)."});
    schema_.add({createdb_opts::MaxSize,
                 "m",
                 "max-size",
                 validate::Required,
                 "  --max-size, -m <size>        \tMaximum database size in GB (default: 500)."});

    schema_.add({0, "", "", option::Arg::None, "\nReporting Options:"});
    add_report_options(schema_);

    schema_.add({0, "", "", option::Arg::None, "\nRuntime Options:"});
    add_runtime_options(schema_, runtime_spec_);

    schema_.add({0, "", "", option::Arg::None, "\nGlobal Options:"});
    add_global_options(schema_);
  }

  [[nodiscard]] std::string_view name() const override { return "createdb"; }

  [[nodiscard]] std::string_view summary() const override { return Summary; }

  [[nodiscard]] const OptionSchema &schema() const override { return schema_; }

  [[nodiscard]] std::any parse_config(const ParsedArgs &args) const override {
    CreateDbConfig config;
    config.source  = parse_source_config(args, source_spec_);
    config.runtime = parse_runtime_config(args, runtime_spec_);
    config.report  = parse_report_config(args);

    if (args.has(createdb_opts::OutputPath)) {
      config.database_path = args.get_string(createdb_opts::OutputPath);
      if (config.database_path.empty()) {
        throw CliUsageError("--output requires a value.");
      }
    } else {
      config.database_path = "createdb_" + current_timestamp_string();
    }

    if (args.has(createdb_opts::MaxSize)) {
      config.max_size_gb = parse_size_option(args.get_string(createdb_opts::MaxSize), "--max-size");
    }

    return std::make_any<CreateDbConfig>(std::move(config));
  }

  [[nodiscard]] PipelinePlan build_plan(const std::any &config) const override {
    CreateDbConfig cfg = std::any_cast<const CreateDbConfig &>(config);
    PipelinePlan plan;
    plan.report_label      = "createdb";
    plan.reporter          = cfg.report.reporter;
    plan.save_run_report   = cfg.report.save_run_report;
    plan.run_report_prefix = "createdb";
    plan.threads           = static_cast<std::size_t>(cfg.runtime.threads);
    plan.success_message   = "Database creation completed successfully!";

    Logger::get_logger()->info("Creating database...");
    switch (cfg.source.mode) {
      case SourceConfig::Mode::Directory:
        Logger::get_logger()->info("Source directory: {}", cfg.source.directory_path);
        Logger::get_logger()->info("Extensions: {}", describe_extensions(cfg.source.extensions));
        Logger::get_logger()->info("Recursive: {}", cfg.source.recursive ? "Yes" : "No");
        break;
      case SourceConfig::Mode::Vector:
        Logger::get_logger()->info("Source files: {} file(s)", cfg.source.file_vector.size());
        break;
      case SourceConfig::Mode::FileList:
        Logger::get_logger()->info("Source file list: {}", cfg.source.file_list_path);
        break;
      case SourceConfig::Mode::Database:
        break;
    }
    Logger::get_logger()->info("Writing to: {}/", cfg.database_path);
    Logger::get_logger()->info("Batch size: {}", cfg.runtime.batch_size);
    Logger::get_logger()->info("Threads: {}", cfg.runtime.threads);
    Logger::get_logger()->info("Max size: {} GB", cfg.max_size_gb);
    plan.output_files.push_back(cfg.database_path + "/");

    auto db = std::make_shared<LMDBDatabase>(cfg.database_path, cfg.max_size_gb);

    plan.source_factory = [cfg]() -> PipelinePlan::SourcePtr {
      using Mode = SourceConfig::Mode;
      switch (cfg.source.mode) {
        case Mode::Directory: {
          auto source = P::from_directory(cfg.source.directory_path,
                                          cfg.source.extensions,
                                          cfg.source.recursive,
                                          cfg.runtime.batch_size);
          return PipelinePlan::SourcePtr(std::move(source));
        }
        case Mode::Vector: {
          auto source = P::from_vector(cfg.source.file_vector);
          return PipelinePlan::SourcePtr(std::move(source));
        }
        case Mode::FileList: {
          auto source = P::from_filelist(cfg.source.file_list_path);
          return PipelinePlan::SourcePtr(std::move(source));
        }
        case Mode::Database:
          break;
      }
      throw CliUsageError("createdb does not support database sources");
    };

    PipelineTask task;
    task.name        = "createdb";
    task.task        = std::make_shared<A::ModelPackTask>("db");
    task.thread_safe = true;
    plan.tasks.push_back(std::move(task));

    PipelineSink sink;
    sink.channel = "db";
    sink.sink    = std::make_shared<P::LmdbSink>(db, cfg.runtime.batch_size);
    plan.sinks.push_back(std::move(sink));

    return plan;
  }

private:
  OptionSchema schema_;
  SourceOptionSpec source_spec_;
  RuntimeOptionSpec runtime_spec_;
};

} // namespace

const CommandSpec &get_createdb_spec() noexcept {
  static const CreateDbSpec spec;
  return spec;
}

} // namespace lahuta::cli
