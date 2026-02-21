/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto t1 = std::make_tuple("besian");
 *   auto t2 = std::make_tuple("sejdiu");
 *   auto t3 = std::make_tuple("@gmail.com");
 *   auto combined = std::tuple_cat(t1, t2, t3);
 *   return std::apply([](auto... args) {
 *     return (std::string{} + ... + std::string(args));
 *   }, combined);
 * }();
 *
 */

#include <cmath>
#include <memory>
#include <string>
#include <string_view>
#include <utility>

#include "analysis/system/model_parse_task.hpp"
#include "analysis/sasa/computation.hpp"
#include "logging/logging.hpp"
#include "parsing/arg_validation.hpp"
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

constexpr std::string_view Summary = "Compute Shrake-Rupley solvent accessible surface area.";

namespace sasa_sr_opts {
constexpr unsigned BaseIndex = 200;
enum : unsigned { //
  Output = BaseIndex,
  OutputStdout,
  ProbeRadius,
  Points,
  IncludeTotal,
  ShowAtomInfo,
  UseBitmask,
  NoSimd
};
} // namespace sasa_sr_opts

struct SasaSrCliConfig {
  SourceConfig source;
  RuntimeConfig runtime;
  ReportConfig report;
  bool output_stdout    = false;
  std::string output_path;
  double probe_radius   = 1.4;
  std::size_t n_points  = 128;
  bool include_total    = false;
  bool show_atom_info   = false;
  bool use_bitmask      = false;
  bool use_simd         = true;
  P::SasaSrParams params;
};

class SasaSrSpec final : public CommandSpec {
public:
  SasaSrSpec() {
    source_spec_.include_is_af2_model = true;
    source_spec_.default_extensions   = {".cif", ".cif.gz", ".pdb", ".pdb.gz"};

    runtime_spec_.include_writer_threads = true;
    runtime_spec_.default_batch_size     = 512;
    runtime_spec_.default_writer_threads = P::get_default_backpressure_config().writer_threads;

    schema_.add({0,
                 "",
                 "",
                 validate::Unknown,
                 std::string("Usage: lahuta sasa-sr [--output <file>] [options]\n"
                             "Author: ")
                     .append(Author)
                     .append("\n\n")
                     .append(Summary)
                     .append("\n"
                             "Outputs: SASA values in JSONL format.\n"
                             "Default output: {\"model\":\"...\",\"sasa\":[...]} with per-atom SASA values.\n"
                             "Use --show-atom-info for labeled output with atom identifiers.\n"
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
    schema_.add({sasa_sr_opts::Output,
                 "o",
                 "output",
                 validate::Required,
                 "  --output, -o <file>          \tOutput file for SASA JSONL (default: sasa_sr.jsonl). "
                 "Use '-' for stdout."});
    schema_.add({sasa_sr_opts::OutputStdout,
                 "",
                 "stdout",
                 option::Arg::None,
                 "  --stdout                     \tWrite JSONL to stdout (same as --output -)."});
    schema_.add({sasa_sr_opts::IncludeTotal,
                 "",
                 "include-total",
                 option::Arg::None,
                 "  --include-total              \tInclude total SASA in JSON output."});
    schema_.add({sasa_sr_opts::ShowAtomInfo,
                 "",
                 "show-atom-info",
                 option::Arg::None,
                 "  --show-atom-info             \tInclude atom labels in output (default: values only)."});

    schema_.add({0, "", "", option::Arg::None, "\nCompute Options:"});
    schema_.add({sasa_sr_opts::ProbeRadius,
                 "",
                 "probe-radius",
                 validate::Required,
                 "  --probe-radius <r>           \tProbe radius in angstroms (default: 1.4)."});
    schema_.add({sasa_sr_opts::Points,
                 "",
                 "points",
                 validate::Required,
                 "  --points <N>                 \tNumber of sphere points per atom (default: 128; "
                 "bitmask fast path uses 64/128/256)."});
    schema_.add({sasa_sr_opts::UseBitmask,
                 "",
                 "use-bitmask",
                 option::Arg::None,
                 "  --use-bitmask                \tEnable bitmask fast path (64/128/256 points)."});
    schema_.add({sasa_sr_opts::NoSimd,
                 "",
                 "no-simd",
                 option::Arg::None,
                 "  --no-simd                    \tDisable SIMD optimizations for bitmask."});

    schema_.add({0, "", "", option::Arg::None, "\nReporting Options:"});
    add_report_options(schema_);

    schema_.add({0, "", "", option::Arg::None, "\nRuntime Options:"});
    add_runtime_options(schema_, runtime_spec_);

    schema_.add({0, "", "", option::Arg::None, "\nGlobal Options:"});
    add_global_options(schema_);
  }

  [[nodiscard]] std::string_view name() const override { return "sasa-sr"; }

  [[nodiscard]] std::string_view summary() const override { return Summary; }

  [[nodiscard]] const OptionSchema &schema() const override { return schema_; }

  [[nodiscard]] std::any parse_config(const ParsedArgs &args) const override {
    SasaSrCliConfig config;
    config.source  = parse_source_config(args, source_spec_);
    config.runtime = parse_runtime_config(args, runtime_spec_);
    config.report  = parse_report_config(args);

    if (config.source.mode == SourceConfig::Mode::Database) {
      config.source.is_af2_model = true;
    }

    if (config.source.mode != SourceConfig::Mode::Database && !config.source.is_af2_model) {
      throw CliUsageError("sasa-sr expects AlphaFold2 model inputs. For file-based sources, pass "
                          "--is_af2_model (or use --database).");
    }

    config.output_stdout = args.get_flag(sasa_sr_opts::OutputStdout);

    if (args.has(sasa_sr_opts::Output)) {
      const auto output_arg = args.get_string(sasa_sr_opts::Output);
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
      config.output_path = "sasa_sr.jsonl";
    }

    if (args.has(sasa_sr_opts::ProbeRadius)) {
      config.probe_radius = std::stod(args.get_string(sasa_sr_opts::ProbeRadius));
      if (!std::isfinite(config.probe_radius) || config.probe_radius < 0.0) {
        throw CliUsageError("--probe-radius must be finite and >= 0.");
      }
    }

    if (args.has(sasa_sr_opts::Points)) {
      const auto raw = std::stoll(args.get_string(sasa_sr_opts::Points));
      if (raw <= 0) {
        throw CliUsageError("--points must be > 0.");
      }
      config.n_points = static_cast<std::size_t>(raw);
    }

    config.include_total  = args.get_flag(sasa_sr_opts::IncludeTotal);
    config.show_atom_info = args.get_flag(sasa_sr_opts::ShowAtomInfo);
    config.use_bitmask    = args.get_flag(sasa_sr_opts::UseBitmask);
    config.use_simd       = !args.get_flag(sasa_sr_opts::NoSimd);

    P::SasaSrParams params;
    params.params.probe_radius = config.probe_radius;
    params.params.n_points     = config.n_points;
    params.params.use_bitmask  = config.use_bitmask;
    params.params.use_simd     = config.use_simd;
    params.include_total       = config.include_total;
    params.show_atom_info      = config.show_atom_info;
    params.channel             = std::string(A::SasaSrOutputChannel);
    params.counters            = std::make_shared<A::SasaSrCounters>();
    config.params              = std::move(params);

    return std::make_any<SasaSrCliConfig>(std::move(config));
  }

  [[nodiscard]] PipelinePlan build_plan(const std::any &config) const override {
    const auto &cfg = std::any_cast<const SasaSrCliConfig &>(config);
    PipelinePlan plan;
    plan.report_label      = "sasa-sr";
    plan.reporter          = cfg.report.reporter;
    plan.save_run_report   = cfg.report.save_run_report;
    plan.run_report_prefix = "sasa-sr";
    plan.threads           = static_cast<std::size_t>(cfg.runtime.threads);
    plan.success_message   = "SASA-SR computation completed successfully!";

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
      throw CliUsageError("sasa-sr does not support this source mode");
    };

    const bool needs_parse_task = cfg.source.mode != SourceConfig::Mode::Database;
    if (needs_parse_task) {
      PipelineTask parse_task;
      parse_task.name        = "parse_model";
      parse_task.task        = std::make_shared<A::ModelParseTask>();
      parse_task.thread_safe = true;
      plan.tasks.push_back(std::move(parse_task));
    }

    PipelineComputation sasa_task;
    sasa_task.name = "sasa_sr";
    if (needs_parse_task) {
      sasa_task.deps.emplace_back("parse_model");
    }
    sasa_task.factory = [params = cfg.params]() {
      return std::make_unique<A::SasaSrComputation>("sasa_sr", params);
    };
    sasa_task.thread_safe = true;
    plan.computations.push_back(std::move(sasa_task));

    auto sink_cfg           = P::get_default_backpressure_config();
    sink_cfg.writer_threads = cfg.runtime.writer_threads;

    PipelineSink data_sink;
    data_sink.channel      = cfg.params.channel;
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

    Logger::get_logger()->info("SASA-SR probe radius: {}", cfg.params.params.probe_radius);
    Logger::get_logger()->info("SASA-SR points: {}", cfg.params.params.n_points);
    return plan;
  }

private:
  OptionSchema schema_;
  SourceOptionSpec source_spec_;
  RuntimeOptionSpec runtime_spec_;
};

} // namespace

const CommandSpec &get_sasa_sr_spec() noexcept {
  static const SasaSrSpec spec;
  return spec;
}

} // namespace lahuta::cli
