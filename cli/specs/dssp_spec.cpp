#include <memory>
#include <string>
#include <string_view>
#include <utility>

#include "analysis/dssp/computation.hpp"
#include "analysis/dssp/records.hpp"
#include "analysis/system/model_parse_task.hpp"
#include "logging/logging.hpp"
#include "parsing/arg_validation.hpp"
#include "parsing/extension_utils.hpp"
#include "parsing/usage_error.hpp"
#include "pipeline/ingest/factory.hpp"
#include "pipeline/runtime/api.hpp"
#include "schemas/shared_options.hpp"
#include "sinks/logging.hpp"
#include "sinks/ndjson.hpp"
#include "specs/command_spec.hpp"

namespace lahuta::cli {
namespace A = lahuta::analysis;
namespace P = lahuta::pipeline;
namespace {

constexpr std::string_view Summary = "Compute DSSP secondary-structure assignments.";

namespace dssp_opts {
constexpr unsigned BaseIndex = 240;
enum : unsigned { //
  Output = BaseIndex,
  OutputStdout,
  NoPreferPi,
  PpStretchLength
};
} // namespace dssp_opts

struct DsspCliConfig {
  SourceConfig source;
  RuntimeConfig runtime;
  ReportConfig report;
  bool output_stdout = false;
  std::string output_path;
  P::DsspParams params;
};

class DsspSpec final : public CommandSpec {
public:
  DsspSpec() {
    source_spec_.include_is_af2_model = true;
    source_spec_.default_extensions   = {".cif", ".cif.gz", ".pdb", ".pdb.gz"};

    runtime_spec_.include_writer_threads = true;
    runtime_spec_.default_batch_size     = 512;
    runtime_spec_.default_writer_threads = P::get_default_backpressure_config().writer_threads;

    schema_.add({0,
                 "",
                 "",
                 validate::Unknown,
                 std::string("Usage: lahuta dssp [--output <file>] [options]\n"
                             "Author: ")
                     .append(Author)
                     .append("\n\n")
                     .append(Summary)
                     .append("\n"
                             "Outputs: DSSP assignments in JSONL format.\n")});

    schema_.add({0, "", "", option::Arg::None, "\nInput Options (choose one):"});
    add_source_options(schema_, source_spec_);

    schema_.add({0, "", "", option::Arg::None, "\nOutput Options:"});
    schema_.add({dssp_opts::Output,
                 "o",
                 "output",
                 validate::Required,
                 "  --output, -o <file>          \tOutput file for DSSP JSONL (default: dssp.jsonl). "
                 "Use '-' for stdout."});
    schema_.add({dssp_opts::OutputStdout,
                 "",
                 "stdout",
                 option::Arg::None,
                 "  --stdout                     \tWrite JSONL to stdout (same as --output -)."});

    schema_.add({0, "", "", option::Arg::None, "\nCompute Options:"});
    schema_.add({dssp_opts::NoPreferPi,
                 "",
                 "no-prefer-pi",
                 option::Arg::None,
                 "  --no-prefer-pi               \tDo not allow pi helices to override alpha helices."});
    schema_.add({dssp_opts::PpStretchLength,
                 "",
                 "pp-stretch-length",
                 validate::Required,
                 "  --pp-stretch-length <n>      \tPPII helix stretch length (2 or 3, default: 2)."});

    schema_.add({0, "", "", option::Arg::None, "\nReporting Options:"});
    add_report_options(schema_);

    schema_.add({0, "", "", option::Arg::None, "\nRuntime Options:"});
    add_runtime_options(schema_, runtime_spec_);

    schema_.add({0, "", "", option::Arg::None, "\nGlobal Options:"});
    add_global_options(schema_);
  }

  [[nodiscard]] std::string_view name() const override { return "dssp"; }
  [[nodiscard]] std::string_view summary() const override { return Summary; }
  [[nodiscard]] const OptionSchema &schema() const override { return schema_; }

  [[nodiscard]] std::any parse_config(const ParsedArgs &args) const override {
    DsspCliConfig config;
    config.source  = parse_source_config(args, source_spec_);
    config.runtime = parse_runtime_config(args, runtime_spec_);
    config.report  = parse_report_config(args);

    if (config.source.mode == SourceConfig::Mode::Database) {
      config.source.is_af2_model = true;
    }

    config.output_stdout = args.get_flag(dssp_opts::OutputStdout);

    if (args.has(dssp_opts::Output)) {
      const auto output_arg = args.get_string(dssp_opts::Output);
      if (output_arg.empty()) throw CliUsageError("--output requires a value.");
      if (output_arg == "-") {
        config.output_stdout = true;
      } else {
        config.output_path = ensure_jsonl_extension(output_arg);
      }
    }

    if (config.output_stdout && !config.output_path.empty()) {
      throw CliUsageError("--stdout cannot be combined with --output <path>. Use --output - for stdout.");
    }

    if (!config.output_stdout && config.output_path.empty()) {
      config.output_path = "dssp.jsonl";
    }

    P::DsspParams params;
    params.channel           = std::string(A::DsspOutputChannel);
    params.prefer_pi_helices = !args.get_flag(dssp_opts::NoPreferPi);

    if (args.has(dssp_opts::PpStretchLength)) {
      const auto raw = std::stoll(args.get_string(dssp_opts::PpStretchLength));
      if (raw != 2 && raw != 3) {
        throw CliUsageError("--pp-stretch-length must be 2 or 3.");
      }
      params.pp_stretch_length = static_cast<int>(raw);
    }

    config.params = std::move(params);

    return std::make_any<DsspCliConfig>(std::move(config));
  }

  [[nodiscard]] PipelinePlan build_plan(const std::any &config) const override {
    const auto &cfg = std::any_cast<const DsspCliConfig &>(config);
    PipelinePlan plan;
    plan.report_label      = "dssp";
    plan.reporter          = cfg.report.reporter;
    plan.save_run_report   = cfg.report.save_run_report;
    plan.run_report_prefix = "dssp";
    plan.threads           = static_cast<std::size_t>(cfg.runtime.threads);
    plan.success_message   = "DSSP computation completed successfully!";

    const bool use_model_pipeline = (cfg.source.is_af2_model ||
                                     cfg.source.mode == SourceConfig::Mode::Database);

    plan.auto_builtins = !use_model_pipeline;
    if (!use_model_pipeline) {
      plan.override_system_params = true;
      plan.system_params.is_model = false;
    }

    if (cfg.source.mode == SourceConfig::Mode::Directory) {
      Logger::get_logger()->info("Source directory: {}", cfg.source.directory_path);
      Logger::get_logger()->info("Extensions: {}", describe_extensions(cfg.source.extensions));
      Logger::get_logger()->info("Recursive: {}", cfg.source.recursive ? "Yes" : "No");
    }

    plan.source_factory = [cfg]() -> PipelinePlan::SourcePtr {
      switch (cfg.source.mode) {
        case SourceConfig::Mode::Database: {
          auto source = P::from_lmdb(cfg.source.database_path,
                                     std::string{},
                                     cfg.runtime.batch_size,
                                     {static_cast<std::size_t>(cfg.runtime.threads) + 1});
          return PipelinePlan::SourcePtr(std::move(source));
        }
        case SourceConfig::Mode::Directory: {
          auto source = P::from_directory(cfg.source.directory_path,
                                          cfg.source.extensions,
                                          cfg.source.recursive,
                                          cfg.runtime.batch_size);
          return PipelinePlan::SourcePtr(std::move(source));
        }
        case SourceConfig::Mode::Vector: {
          auto source = P::from_vector(cfg.source.file_vector);
          return PipelinePlan::SourcePtr(std::move(source));
        }
        case SourceConfig::Mode::FileList: {
          auto source = P::from_filelist(cfg.source.file_list_path);
          return PipelinePlan::SourcePtr(std::move(source));
        }
      }
      throw CliUsageError("dssp does not support this source mode");
    };

    const bool needs_parse_task = use_model_pipeline && cfg.source.mode != SourceConfig::Mode::Database;
    if (needs_parse_task) {
      PipelineTask parse_task;
      parse_task.name        = "parse_model";
      parse_task.task        = std::make_shared<A::ModelParseTask>();
      parse_task.thread_safe = true;
      plan.tasks.push_back(std::move(parse_task));
    }

    PipelineComputation computation;
    computation.name = "dssp";
    if (use_model_pipeline) {
      computation.factory = [params = cfg.params]() {
        return std::make_unique<A::DsspModelComputation>("dssp", params);
      };
    } else {
      computation.factory = [params = cfg.params]() {
        return std::make_unique<A::DsspComputation>("dssp", params);
      };
    }
    if (needs_parse_task) {
      computation.deps.emplace_back("parse_model");
    }
    computation.thread_safe = true;
    plan.computations.push_back(std::move(computation));

    auto sink_cfg           = P::get_default_backpressure_config();
    sink_cfg.writer_threads = cfg.runtime.writer_threads;

    PipelineSink sink;
    sink.channel      = cfg.params.channel;
    sink.backpressure = sink_cfg;
    if (cfg.output_stdout) {
      sink.sink = std::make_shared<P::LoggingSink>();
      Logger::get_logger()->info("Writing to: stdout");
    } else {
      sink.sink = std::make_shared<P::NdjsonFileSink>(cfg.output_path);
      Logger::get_logger()->info("Writing to: {}", cfg.output_path);
      plan.output_files.push_back(cfg.output_path);
    }
    plan.sinks.push_back(std::move(sink));

    return plan;
  }

private:
  OptionSchema schema_;
  SourceOptionSpec source_spec_;
  RuntimeOptionSpec runtime_spec_;
};

} // namespace

const CommandSpec &get_dssp_spec() noexcept {
  static const DsspSpec spec;
  return spec;
}

} // namespace lahuta::cli
