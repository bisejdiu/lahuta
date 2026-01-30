#include <cmath>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "analysis/extract/extract_tasks.hpp"
#include "analysis/sasa/computation.hpp"
#include "logging/logging.hpp"
#include "parsing/arg_validation.hpp"
#include "pipeline/ingest/factory.hpp"
#include "pipeline/runtime/api.hpp"
#include "schemas/shared_options.hpp"
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
  OutputDir = BaseIndex,
  ProbeRadius,
  Points,
  IncludeTotal,
  UseBitmask,
  NoSimd
};
} // namespace sasa_sr_opts

struct SasaSrCliConfig {
  SourceConfig source;
  RuntimeConfig runtime;
  ReportConfig report;
  std::filesystem::path output_dir;
  double probe_radius  = 1.4;
  std::size_t n_points = 128;
  bool include_total   = false;
  bool use_bitmask     = false;
  bool use_simd        = true;
  P::SasaSrParams params;
};

class SasaSrSpec final : public CommandSpec {
public:
  SasaSrSpec() {
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
                 std::string("Usage: lahuta sasa-sr [--output-dir <dir>] [options]\n\n")
                     .append(Summary)
                     .append("\n"
                             "Outputs: per_protein_sasa_sr.jsonl (NDJSON) in the output directory.\n"
                             "Each entry has an \"Atom\" array of {label: sasa} objects where label is\n"
                             "atom-index-atom_name-residue_id-residue_name-chain_id.\n"
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
    schema_.add({sasa_sr_opts::OutputDir,
                 "o",
                 "output-dir",
                 validate::Required,
                 "  --output-dir, -o <dir>       \tOutput directory for per-model JSON files (default: .)."});
    schema_.add({sasa_sr_opts::IncludeTotal,
                 "",
                 "include-total",
                 option::Arg::None,
                 "  --include-total              \tInclude total SASA in JSON output."});

    schema_.add({0, "", "", option::Arg::None, "\nAlgorithm Options:"});
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
      throw std::runtime_error("sasa-sr expects AlphaFold2 model inputs. For file-based sources, pass "
                               "--is_af2_model (or use --database).");
    }

    std::string output_arg;
    if (args.has(sasa_sr_opts::OutputDir)) {
      output_arg = args.get_string(sasa_sr_opts::OutputDir);
      if (output_arg.empty()) {
        throw std::runtime_error("--output-dir requires a value.");
      }
      config.output_dir = std::filesystem::path(output_arg);
    } else {
      config.output_dir = std::filesystem::path(".");
      output_arg        = config.output_dir.string();
    }

    if (args.has(sasa_sr_opts::ProbeRadius)) {
      config.probe_radius = std::stod(args.get_string(sasa_sr_opts::ProbeRadius));
      if (!std::isfinite(config.probe_radius) || config.probe_radius < 0.0) {
        throw std::runtime_error("--probe-radius must be finite and >= 0.");
      }
    }

    if (args.has(sasa_sr_opts::Points)) {
      const auto raw = std::stoll(args.get_string(sasa_sr_opts::Points));
      if (raw <= 0) {
        throw std::runtime_error("--points must be > 0.");
      }
      config.n_points = static_cast<std::size_t>(raw);
    }

    config.include_total = args.get_flag(sasa_sr_opts::IncludeTotal);
    config.use_bitmask   = args.get_flag(sasa_sr_opts::UseBitmask);
    config.use_simd      = !args.get_flag(sasa_sr_opts::NoSimd);
    if (config.output_dir.empty()) {
      throw std::runtime_error("Output directory cannot be empty.");
    }

    std::error_code ec;
    if (std::filesystem::exists(config.output_dir, ec)) {
      if (ec) {
        throw std::runtime_error("Unable to access output directory: " + output_arg);
      }
      if (!std::filesystem::is_directory(config.output_dir, ec)) {
        throw std::runtime_error("Output path is not a directory: " + output_arg);
      }
    } else {
      if (!std::filesystem::create_directories(config.output_dir, ec) || ec) {
        throw std::runtime_error("Unable to create output directory: " + output_arg);
      }
    }

    P::SasaSrParams params;
    params.params.probe_radius = config.probe_radius;
    params.params.n_points     = config.n_points;
    params.params.use_bitmask  = config.use_bitmask;
    params.params.use_simd     = config.use_simd;
    params.include_total       = config.include_total;
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
      throw std::runtime_error("sasa-sr does not support this source mode");
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
    const auto output_path = (cfg.output_dir / "per_protein_sasa_sr.jsonl").string();
    data_sink.sink         = std::make_shared<P::NdjsonFileSink>(output_path);
    Logger::get_logger()->info("SASA-SR output -> {}", output_path);
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
