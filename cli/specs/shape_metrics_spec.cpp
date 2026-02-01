#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <utility>

#include "analysis/extract/extract_tasks.hpp"
#include "logging/logging.hpp"
#include "parsing/arg_validation.hpp"
#include "parsing/usage_error.hpp"
#include "pipeline/runtime/api.hpp"
#include "schemas/shared_options.hpp"
#include "sinks/ndjson.hpp"
#include "specs/command_spec.hpp"
#include "tasks/shape_metrics_summary_sink.hpp"
#include "tasks/shape_metrics_task.hpp"

namespace lahuta::cli {
namespace A = lahuta::analysis;
namespace P = lahuta::pipeline;
namespace {

constexpr std::string_view Summary = "Compute tensor-based shape metrics with pLDDT/DSSP trimming.";

namespace shape_metrics_opts {
constexpr unsigned BaseIndex = 200;
enum : unsigned { OutputDir = BaseIndex, MinHighFraction };
} // namespace shape_metrics_opts

struct ShapeMetricsCliConfig {
  SourceConfig source;
  RuntimeConfig runtime;
  ReportConfig report;
  std::filesystem::path output_dir;
  double min_high_fraction = 0.80;
  std::shared_ptr<const shape_metrics::ShapeMetricsConfig> task_config;
};

class ShapeMetricsSpec final : public CommandSpec {
public:
  ShapeMetricsSpec() {
    source_spec_.include_is_af2_model = true;
    source_spec_.default_extensions   = {".cif", ".cif.gz", ".pdb", ".pdb.gz"};

    runtime_spec_.include_writer_threads = true;
    runtime_spec_.default_threads        = 8;
    runtime_spec_.default_batch_size     = 512;
    runtime_spec_.default_writer_threads = P::get_default_backpressure_config().writer_threads;

    schema_.add({0,
                 "",
                 "",
                 validate::Unknown,
                 std::string("Usage: lahuta shape-metrics [--output-dir <dir>] [options]\n"
                             "Author: ")
                     .append(Author)
                     .append("\n\n")
                     .append(Summary)
                     .append("\n"
                             "Metrics include Rg, asphericity, and kappa^2. Trimming removes\n"
                             "low-confidence coil/turn/bend tails.\n\n"
                             "Outputs: per_protein_shape_metrics.jsonl (NDJSON) in the output directory.\n"
                             "Note: file-based inputs must be AF2 model files (mmCIF).")});

    schema_.add({0, "", "", option::Arg::None, "\nInput Options (choose one):"});
    schema_.add({shared_opts::SourceDatabase,
                 "",
                 "database",
                 validate::Required,
                 "  --database <path>            \tProcess structures from database."});
    schema_.add({shared_opts::SourceDirectory,
                 "d",
                 "directory",
                 validate::Required,
                 "  --directory, -d <path>       \tProcess all files in directory."});
    schema_.add({shared_opts::SourceVector,
                 "f",
                 "files",
                 validate::Required,
                 "  --files, -f <file1,file2>    \tProcess specific files (comma-separated or repeat -f)."});
    schema_.add({shared_opts::SourceFileList,
                 "l",
                 "file-list",
                 validate::Required,
                 "  --file-list, -l <path>       \tProcess files listed in text file (one per line)."});

    schema_.add({0, "", "", option::Arg::None, "\nDirectory Options:"});
    schema_.add({shared_opts::SourceExtension,
                 "e",
                 "extension",
                 validate::Required,
                 "  --extension, -e <ext>        \tFile extension(s) for directory mode. Repeat or "
                 "comma-separate values (default: .cif, .cif.gz, .pdb, .pdb.gz)."});
    schema_.add({shared_opts::SourceRecursive,
                 "r",
                 "recursive",
                 option::Arg::None,
                 "  --recursive, -r              \tRecursively search subdirectories."});

    schema_.add({0, "", "", option::Arg::None, "\nModel Options:"});
    schema_.add({shared_opts::SourceIsAf2Model,
                 "",
                 "is_af2_model",
                 option::Arg::None,
                 "  --is_af2_model               \tRequired for file-based sources; inputs are AlphaFold2 "
                 "models (AF2-like mmCIF)."});

    schema_.add({0, "", "", option::Arg::None, "\nOutput Options:"});
    schema_.add({shape_metrics_opts::OutputDir,
                 "o",
                 "output-dir",
                 validate::Required,
                 "  --output-dir, -o <dir>       \tOutput directory for NDJSON results (default: .)."});

    schema_.add({0, "", "", option::Arg::None, "\nAlgorithm Options:"});
    schema_.add({shape_metrics_opts::MinHighFraction,
                 "",
                 "min-high-fraction",
                 validate::Required,
                 "  --min-high-fraction <f>      \tMinimum high-confidence fraction after trimming "
                 "(default: 0.80)."});
    add_report_options(schema_);

    schema_.add({0, "", "", option::Arg::None, "\nRuntime Options:"});
    add_runtime_options(schema_, runtime_spec_);

    schema_.add({0, "", "", option::Arg::None, "\nGlobal Options:"});
    add_global_options(schema_);
  }

  [[nodiscard]] std::string_view name() const override { return "shape-metrics"; }

  [[nodiscard]] std::string_view summary() const override { return Summary; }

  [[nodiscard]] const OptionSchema &schema() const override { return schema_; }

  [[nodiscard]] std::any parse_config(const ParsedArgs &args) const override {
    ShapeMetricsCliConfig config;
    config.source  = parse_source_config(args, source_spec_);
    config.runtime = parse_runtime_config(args, runtime_spec_);
    config.report  = parse_report_config(args);

    if (config.source.mode == SourceConfig::Mode::Database) {
      config.source.is_af2_model = true;
    }

    if (config.source.mode != SourceConfig::Mode::Database && !config.source.is_af2_model) {
      throw CliUsageError("shape-metrics expects AlphaFold2 model inputs. For file-based sources, pass "
                          "--is_af2_model (or use --database).");
    }

    std::string output_arg;
    if (args.has(shape_metrics_opts::OutputDir)) {
      output_arg = args.get_string(shape_metrics_opts::OutputDir);
      if (output_arg.empty()) {
        throw CliUsageError("--output-dir requires a value.");
      }
    } else {
      output_arg = ".";
    }

    if (args.has(shape_metrics_opts::MinHighFraction)) {
      config.min_high_fraction = std::stod(args.get_string(shape_metrics_opts::MinHighFraction));
      if (config.min_high_fraction < 0.0 || config.min_high_fraction > 1.0) {
        throw CliUsageError("--min-high-fraction must be between 0 and 1.");
      }
    }

    config.output_dir = validate_output_dir(output_arg);

    auto counters               = std::make_shared<shape_metrics::ShapeMetricsCounters>();
    auto task_cfg               = std::make_shared<shape_metrics::ShapeMetricsConfig>();
    task_cfg->min_high_fraction = config.min_high_fraction;
    task_cfg->counters          = std::move(counters);
    config.task_config          = std::move(task_cfg);

    return std::make_any<ShapeMetricsCliConfig>(std::move(config));
  }

  [[nodiscard]] PipelinePlan build_plan(const std::any &config) const override {
    const auto &cfg = std::any_cast<const ShapeMetricsCliConfig &>(config);
    PipelinePlan plan;
    plan.report_label      = "shape-metrics";
    plan.reporter          = cfg.report.reporter;
    plan.save_run_report   = cfg.report.save_run_report;
    plan.run_report_prefix = "shape-metrics";
    plan.threads           = static_cast<std::size_t>(cfg.runtime.threads);

    plan.source_factory = [cfg]() -> PipelinePlan::SourcePtr {
      using Mode = SourceConfig::Mode;
      switch (cfg.source.mode) {
        case Mode::Database: {
          auto source = P::from_lmdb(cfg.source.database_path,
                                     std::string{},
                                     cfg.runtime.batch_size,
                                     {static_cast<std::size_t>(cfg.runtime.threads) + 1});
          return PipelinePlan::SourcePtr(std::move(source));
        }
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
      }
      throw CliUsageError("shape-metrics does not support this source mode");
    };

    auto sink_cfg           = P::get_default_backpressure_config();
    sink_cfg.writer_threads = cfg.runtime.writer_threads;

    const bool needs_parse_task = cfg.source.mode != SourceConfig::Mode::Database;
    if (needs_parse_task) {
      PipelineTask parse_task;
      parse_task.name        = "parse_model";
      parse_task.task        = std::make_shared<A::ModelParseTask>();
      parse_task.thread_safe = true;
      plan.tasks.push_back(std::move(parse_task));
    }

    PipelineTask metrics_task;
    metrics_task.name = "shape_metrics";
    metrics_task.task = shape_metrics::make_shape_metrics_task(cfg.task_config);
    if (needs_parse_task) {
      metrics_task.deps.emplace_back("parse_model");
    }
    metrics_task.thread_safe = true;
    plan.tasks.push_back(std::move(metrics_task));

    PipelineSink data_sink;
    data_sink.channel      = std::string(shape_metrics::OutputChannel);
    data_sink.backpressure = sink_cfg;
    const auto output_path = (cfg.output_dir / "per_protein_shape_metrics.jsonl").string();
    data_sink.sink         = std::make_shared<P::NdjsonFileSink>(output_path);
    Logger::get_logger()->info("Shape-metrics output -> {}", output_path);
    plan.sinks.push_back(std::move(data_sink));

    PipelineSink summary_sink;
    summary_sink.channel      = std::string(shape_metrics::OutputChannel);
    summary_sink.backpressure = sink_cfg;
    summary_sink.sink         = shape_metrics::make_shape_metrics_summary_sink(cfg.output_dir,
                                                                       cfg.task_config->counters);
    Logger::get_logger()->info("Shape-metrics summary -> {}",
                               (cfg.output_dir / "shape_metrics_summary.json").string());
    plan.sinks.push_back(std::move(summary_sink));

    return plan;
  }

private:
  OptionSchema schema_;
  SourceOptionSpec source_spec_;
  RuntimeOptionSpec runtime_spec_;
};

} // namespace

const CommandSpec &get_shape_metrics_spec() noexcept {
  static const ShapeMetricsSpec spec;
  return spec;
}

} // namespace lahuta::cli
