#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <utility>

#include "analysis/system/model_parse_task.hpp"
#include "logging/logging.hpp"
#include "parsing/arg_validation.hpp"
#include "parsing/usage_error.hpp"
#include "pipeline/runtime/api.hpp"
#include "schemas/shared_options.hpp"
#include "sinks/logging.hpp"
#include "sinks/ndjson.hpp"
#include "specs/command_spec.hpp"
#include "tasks/compaction_rg_summary_sink.hpp"
#include "tasks/compaction_rg_task.hpp"

namespace lahuta::cli {
namespace A = lahuta::analysis;
namespace P = lahuta::pipeline;
namespace {

constexpr std::string_view Summary = "Compute radius of gyration with pLDDT/DSSP trimming.";

namespace compaction_rg_opts {
constexpr unsigned BaseIndex = 200;
enum : unsigned { Output = BaseIndex, OutputStdout, MinHighFraction };
} // namespace compaction_rg_opts

struct CompactionRgCliConfig {
  SourceConfig source;
  RuntimeConfig runtime;
  ReportConfig report;
  bool output_stdout = false;
  std::string output_path;
  std::filesystem::path summary_path;
  double min_high_fraction = 0.80;
  std::shared_ptr<const compaction_rg::CompactionRgConfig> task_config;
};

class CompactionRgSpec final : public CommandSpec {
public:
  CompactionRgSpec() {
    source_spec_.include_is_af2_model = true;
    source_spec_.default_extensions   = {".cif", ".cif.gz", ".pdb", ".pdb.gz"};

    runtime_spec_.include_writer_threads = true;
    runtime_spec_.default_batch_size     = 512;
    runtime_spec_.default_writer_threads = P::get_default_backpressure_config().writer_threads;

    schema_.add({0,
                 "",
                 "",
                 validate::Unknown,
                 std::string("Usage: lahuta compaction-rg [--output <file>] [options]\n"
                             "Author: ")
                     .append(Author)
                     .append("\n\n")
                     .append(Summary)
                     .append("\n"
                             "Trimming removes N/C termini that are both low-confidence (pLDDT) and\n"
                             "unstructured (DSSP coil/turn/bend).\n\n"
                             "Outputs: Rg metrics in JSONL format (+ summary JSON in same directory).\n"
                             "Note: file-based inputs must be AF2 model files.")});

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

    schema_.add({shared_opts::SourceIsAf2Model,
                 "",
                 "is_af2_model",
                 option::Arg::None,
                 "  --is_af2_model               \tRequired for file-based sources; inputs are AlphaFold2 "
                 "models (AF2-like mmCIF)."});

    schema_.add({0, "", "", option::Arg::None, "\nOutput Options:"});
    schema_.add({compaction_rg_opts::Output,
                 "o",
                 "output",
                 validate::Required,
                 "  --output, -o <file>          \tOutput file for Rg JSONL (default: compaction_rg.jsonl). "
                 "Use '-' for stdout."});
    schema_.add({compaction_rg_opts::OutputStdout,
                 "",
                 "stdout",
                 option::Arg::None,
                 "  --stdout                     \tWrite JSONL to stdout (same as --output -)."});

    schema_.add({0, "", "", option::Arg::None, "\nCompute Options:"});
    schema_.add({compaction_rg_opts::MinHighFraction,
                 "",
                 "min-high-fraction",
                 validate::Required,
                 "  --min-high-fraction <f>      \tMinimum high-confidence fraction after trimming "
                 "(default: 0.80)."});
    schema_.add({0, "", "", option::Arg::None, "\nReporting Options:"});
    add_report_options(schema_);

    schema_.add({0, "", "", option::Arg::None, "\nRuntime Options:"});
    add_runtime_options(schema_, runtime_spec_);

    schema_.add({0, "", "", option::Arg::None, "\nGlobal Options:"});
    add_global_options(schema_);
  }

  [[nodiscard]] std::string_view name() const override { return "compaction-rg"; }

  [[nodiscard]] std::string_view summary() const override { return Summary; }

  [[nodiscard]] const OptionSchema &schema() const override { return schema_; }

  [[nodiscard]] std::any parse_config(const ParsedArgs &args) const override {
    CompactionRgCliConfig config;
    config.source  = parse_source_config(args, source_spec_);
    config.runtime = parse_runtime_config(args, runtime_spec_);
    config.report  = parse_report_config(args);

    if (config.source.mode == SourceConfig::Mode::Database) {
      config.source.is_af2_model = true;
    }

    if (config.source.mode != SourceConfig::Mode::Database && !config.source.is_af2_model) {
      throw CliUsageError("compaction-rg expects AlphaFold2 model inputs. For file-based sources, pass "
                          "--is_af2_model (or use --database).");
    }

    config.output_stdout = args.get_flag(compaction_rg_opts::OutputStdout);

    if (args.has(compaction_rg_opts::Output)) {
      const auto output_arg = args.get_string(compaction_rg_opts::Output);
      if (output_arg.empty()) {
        throw CliUsageError("--output requires a value.");
      }
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
      config.output_path = "compaction_rg.jsonl";
    }

    // Get summary path from output path: foo.jsonl -> foo_summary.json
    // For stdout mode, use default summary path
    if (config.output_stdout) {
      config.summary_path = "compaction_rg_summary.json";
    } else {
      std::filesystem::path output_fs(config.output_path);
      config.summary_path = output_fs.parent_path() / (output_fs.stem().string() + "_summary.json");
    }

    if (args.has(compaction_rg_opts::MinHighFraction)) {
      config.min_high_fraction = std::stod(args.get_string(compaction_rg_opts::MinHighFraction));
      if (config.min_high_fraction < 0.0 || config.min_high_fraction > 1.0) {
        throw CliUsageError("--min-high-fraction must be between 0 and 1.");
      }
    }

    auto counters               = std::make_shared<compaction_rg::CompactionRgCounters>();
    auto task_cfg               = std::make_shared<compaction_rg::CompactionRgConfig>();
    task_cfg->min_high_fraction = config.min_high_fraction;
    task_cfg->counters          = std::move(counters);
    config.task_config          = std::move(task_cfg);

    return std::make_any<CompactionRgCliConfig>(std::move(config));
  }

  [[nodiscard]] PipelinePlan build_plan(const std::any &config) const override {
    const auto &cfg = std::any_cast<const CompactionRgCliConfig &>(config);
    PipelinePlan plan;
    plan.report_label      = "compaction-rg";
    plan.reporter          = cfg.report.reporter;
    plan.save_run_report   = cfg.report.save_run_report;
    plan.run_report_prefix = "compaction-rg";
    plan.threads           = static_cast<std::size_t>(cfg.runtime.threads);
    plan.success_message   = "Compaction-rg computation completed successfully!";

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
      throw CliUsageError("compaction-rg does not support this source mode");
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

    PipelineTask rg_task;
    rg_task.name = "compaction_rg";
    rg_task.task = compaction_rg::make_compaction_rg_task(cfg.task_config);
    if (needs_parse_task) {
      rg_task.deps.emplace_back("parse_model");
    }
    rg_task.thread_safe = true;
    plan.tasks.push_back(std::move(rg_task));

    PipelineSink data_sink;
    data_sink.channel      = std::string(compaction_rg::OutputChannel);
    data_sink.backpressure = sink_cfg;
    if (cfg.output_stdout) {
      data_sink.sink = std::make_shared<P::LoggingSink>();
      Logger::get_logger()->info("Writing to: stdout");
    } else {
      data_sink.sink = std::make_shared<P::NdjsonFileSink>(cfg.output_path);
      Logger::get_logger()->info("Writing to: {}", cfg.output_path);
      plan.output_files.push_back(cfg.output_path);
    }
    plan.sinks.push_back(std::move(data_sink));

    PipelineSink summary_sink;
    summary_sink.channel      = std::string(compaction_rg::OutputChannel);
    summary_sink.backpressure = sink_cfg;
    summary_sink.sink         = compaction_rg::make_compaction_rg_summary_sink(cfg.summary_path,
                                                                       cfg.task_config->counters);
    plan.sinks.push_back(std::move(summary_sink));

    Logger::get_logger()->info("Writing to: {}", cfg.summary_path.string());
    plan.output_files.push_back(cfg.summary_path.string());

    return plan;
  }

private:
  OptionSchema schema_;
  SourceOptionSpec source_spec_;
  RuntimeOptionSpec runtime_spec_;
};

} // namespace

const CommandSpec &get_compaction_rg_spec() noexcept {
  static const CompactionRgSpec spec;
  return spec;
}

} // namespace lahuta::cli
