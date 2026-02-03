/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto parts = std::make_tuple("besian", "sejdiu", "@gmail.com");
 *   return std::apply([](auto... p) { std::string s; (s.append(p), ...); return s; }, parts);
 * }();
 *
 */

#include <memory>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include "analysis/system/model_parse_task.hpp"
#include "logging/logging.hpp"
#include "parsing/arg_validation.hpp"
#include "parsing/extension_utils.hpp"
#include "parsing/usage_error.hpp"
#include "pipeline/runtime/api.hpp"
#include "schemas/shared_options.hpp"
#include "sinks/logging.hpp"
#include "sinks/ndjson.hpp"
#include "specs/command_spec.hpp"
#include "tasks/quality_metrics_task.hpp"

namespace lahuta::cli {
namespace A = lahuta::analysis;
namespace P = lahuta::pipeline;
namespace {

constexpr std::string_view Summary = "Compute per-protein group metrics for pLDDT and DSSP signals.";

namespace quality_metrics_opts {
constexpr unsigned BaseIndex = 200;
enum : unsigned {
  PlddtGroup = BaseIndex,
  DsspGroup,
  SegmentGroup,
  SegmentMin,
  NoOverlap,
  Output,
  OutputStdout
};
} // namespace quality_metrics_opts

struct QualityMetricsCliConfig {
  SourceConfig source;
  RuntimeConfig runtime;
  ReportConfig report;
  bool output_stdout      = false;
  std::string output_path = "per_protein_metrics.jsonl";
  std::shared_ptr<const quality_metrics::QualityMetricsConfig> metrics_config;
};

class QualityMetricsSpec final : public CommandSpec {
public:
  QualityMetricsSpec() {
    source_spec_.include_is_af2_model = true;
    source_spec_.default_extensions   = {".cif", ".cif.gz", ".pdb", ".pdb.gz"};

    runtime_spec_.include_writer_threads = true;
    runtime_spec_.default_batch_size     = 512;
    runtime_spec_.default_writer_threads = P::get_default_backpressure_config().writer_threads;

    schema_.add(
        {0,
         "",
         "",
         validate::Unknown,
         std::string("Usage: lahuta quality-metrics [options]\n"
                     "Author: ")
             .append(Author)
             .append("\n\n")
             .append(Summary)
             .append("\n"
                     "Computes per-protein fractions for user-defined pLDDT and DSSP groups, plus\n"
                     "optional segment statistics for pLDDT groups and overlap fractions between\n"
                     "pLDDT and DSSP groups.\n"
                     "Segment statistics measure how many consecutive residues fall in a pLDDT\n"
                     "group, the longest run, and the fraction in long runs.\n"
                     "Overlap fractions report how often a residue is in both a pLDDT group and a\n"
                     "DSSP group (for example, poor-confidence coil).\n\n"
                     "For example, if you group Low and VeryLow as \"poor\", the output tells\n"
                     "you what fraction of each protein is poor-confidence, plus how much of that\n"
                     "poor-confidence region falls in coils or helices.\n\n"
                     "Valid pLDDT values are: VeryHigh, High, Low, VeryLow.\n"
                     "Valid  DSSP values are: Coil, AlphaHelix, Helix3_10, HelixPi,\n"
                     "                        PolyProlineHelix, Strand, Turn, Bend.\n\n"
                     "Examples:\n"
                     "  lahuta quality-metrics -d dir/ --is_af2_model --output metrics.jsonl\n"
                     "  lahuta quality-metrics -d dir/ --is_af2_model --plddt-group poor=Low,VeryLow\n"
                     "    --segment-group poor\n"
                     "  lahuta quality-metrics --directory /path/to/models --is_af2_model \\\n"
                     "    --plddt-group poor=Low,VeryLow --segment-group poor \\\n"
                     "    --dssp-group coil=Coil,Turn,Bend \\\n"
                     "    --dssp-group helix=AlphaHelix,Helix3_10,HelixPi,PolyProlineHelix \\\n"
                     "    --dssp-group strand=Strand \\\n"
                     "    --output metrics.jsonl")});

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
    schema_.add({quality_metrics_opts::Output,
                 "o",
                 "output",
                 validate::Required,
                 "  --output, -o <path>          \tWrite JSONL to file (default: "
                 "per_protein_metrics.jsonl). Use '-' for stdout."});
    schema_.add({quality_metrics_opts::OutputStdout,
                 "",
                 "stdout",
                 option::Arg::None,
                 "  --stdout                     \tWrite JSONL to stdout (same as --output -)."});

    schema_.add({0, "", "", option::Arg::None, "\nCompute Options:"});
    schema_.add({quality_metrics_opts::PlddtGroup,
                 "",
                 "plddt-group",
                 validate::Required,
                 "  --plddt-group <name=values>  \tAdd a pLDDT group (repeatable)."});
    schema_.add({quality_metrics_opts::DsspGroup,
                 "",
                 "dssp-group",
                 validate::Required,
                 "  --dssp-group <name=values>   \tAdd a DSSP group (repeatable)."});
    schema_.add({quality_metrics_opts::SegmentGroup,
                 "",
                 "segment-group",
                 validate::Required,
                 "  --segment-group <name>       \tEnable segment metrics for a pLDDT group (repeatable)."});
    schema_.add(
        {quality_metrics_opts::SegmentMin,
         "",
         "segment-min",
         validate::Required,
         "  --segment-min <N>            \tMinimum residues for long-segment fraction (default: 30)."});
    schema_.add({quality_metrics_opts::NoOverlap,
                 "",
                 "no-overlap",
                 option::Arg::None,
                 "  --no-overlap                 \tDisable pLDDT x DSSP overlap metrics."});

    schema_.add({0, "", "", option::Arg::None, "\nReporting Options:"});
    add_report_options(schema_);

    schema_.add({0, "", "", option::Arg::None, "\nRuntime Options:"});
    add_runtime_options(schema_, runtime_spec_);

    schema_.add({0, "", "", option::Arg::None, "\nGlobal Options:"});
    add_global_options(schema_);
  }

  [[nodiscard]] std::string_view name() const override { return "quality-metrics"; }

  [[nodiscard]] std::string_view summary() const override { return Summary; }

  [[nodiscard]] const OptionSchema &schema() const override { return schema_; }

  [[nodiscard]] std::any parse_config(const ParsedArgs &args) const override {
    QualityMetricsCliConfig config;
    config.source  = parse_source_config(args, source_spec_);
    config.runtime = parse_runtime_config(args, runtime_spec_);
    config.report  = parse_report_config(args);

    std::unordered_set<std::string> plddt_group_names;
    std::unordered_set<std::string> dssp_group_names;
    std::vector<quality_metrics::GroupSpec> plddt_groups;
    std::vector<quality_metrics::GroupSpec> dssp_groups;

    if (args.has(quality_metrics_opts::PlddtGroup)) {
      const auto raw_values = args.get_all_strings(quality_metrics_opts::PlddtGroup);
      for (const auto &raw : raw_values) {
        quality_metrics::add_group_definition(raw, true, plddt_group_names, plddt_groups);
      }
    }

    if (args.has(quality_metrics_opts::DsspGroup)) {
      const auto raw_values = args.get_all_strings(quality_metrics_opts::DsspGroup);
      for (const auto &raw : raw_values) {
        quality_metrics::add_group_definition(raw, false, dssp_group_names, dssp_groups);
      }
    }

    bool used_default_plddt = false;
    if (plddt_groups.empty()) {
      quality_metrics::add_default_plddt_groups(plddt_groups);
      used_default_plddt = true;
    }
    if (dssp_groups.empty()) {
      quality_metrics::add_default_dssp_groups(dssp_groups);
    }

    std::unordered_set<std::string> segment_group_names;
    if (args.has(quality_metrics_opts::SegmentGroup)) {
      const auto raw_values = args.get_all_strings(quality_metrics_opts::SegmentGroup);
      for (const auto &raw : raw_values) {
        std::vector<std::string> tokens;
        detail::split_argument_list(raw, false, tokens);
        for (const auto &token : tokens) {
          const std::string name = quality_metrics::normalize_group_name(token);
          if (name.empty()) {
            throw CliUsageError("Invalid segment group name '" + token + "'.");
          }
          segment_group_names.insert(name);
        }
      }
    }

    if (segment_group_names.empty() && used_default_plddt) {
      segment_group_names.insert("poor");
    }

    if (!segment_group_names.empty()) {
      std::unordered_set<std::string> available;
      for (const auto &group : plddt_groups) {
        available.insert(group.name);
      }
      for (const auto &name : segment_group_names) {
        if (available.count(name) == 0) {
          throw CliUsageError("Segment group '" + name + "' does not match any pLDDT group.");
        }
      }
      for (auto &group : plddt_groups) {
        group.segment = segment_group_names.count(group.name) > 0;
      }
    }

    quality_metrics::finalize_group_keys(plddt_groups, dssp_groups);

    bool compute_overlaps   = !args.get_flag(quality_metrics_opts::NoOverlap);
    std::size_t segment_min = 30;
    if (args.has(quality_metrics_opts::SegmentMin)) {
      segment_min = std::stoull(args.get_string(quality_metrics_opts::SegmentMin));
      if (segment_min == 0) {
        throw CliUsageError("--segment-min must be positive");
      }
    }

    config.output_stdout = args.get_flag(quality_metrics_opts::OutputStdout);

    if (args.has(quality_metrics_opts::Output)) {
      const std::string output_arg = args.get_string(quality_metrics_opts::Output);
      if (output_arg.empty()) {
        throw CliUsageError("--output requires a value.");
      }
      if (output_arg == "-") {
        config.output_stdout = true;
        config.output_path   = "-";
      } else {
        config.output_path = output_arg;
      }
    }

    if (config.output_stdout && args.has(quality_metrics_opts::Output) && config.output_path != "-") {
      throw CliUsageError("--stdout cannot be combined with --output <path>. Use --output - for stdout.");
    }

    if (config.source.mode == SourceConfig::Mode::Database) {
      config.source.is_af2_model = true;
    }

    if (config.source.mode != SourceConfig::Mode::Database && !config.source.is_af2_model) {
      throw CliUsageError("quality-metrics expects AlphaFold2 model inputs. For file-based sources, "
                          "pass --is_af2_model (or use --database).");
    }

    auto metrics_config              = std::make_shared<quality_metrics::QualityMetricsConfig>();
    metrics_config->plddt_groups     = std::move(plddt_groups);
    metrics_config->dssp_groups      = std::move(dssp_groups);
    metrics_config->compute_overlaps = compute_overlaps;
    metrics_config->segment_min      = segment_min;
    if (metrics_config->compute_overlaps && !metrics_config->plddt_groups.empty() &&
        !metrics_config->dssp_groups.empty()) {
      metrics_config->overlaps = quality_metrics::build_overlap_specs(metrics_config->plddt_groups,
                                                                      metrics_config->dssp_groups);
    }
    config.metrics_config = std::move(metrics_config);

    return std::make_any<QualityMetricsCliConfig>(std::move(config));
  }

  [[nodiscard]] PipelinePlan build_plan(const std::any &config) const override {
    const auto &cfg = std::any_cast<const QualityMetricsCliConfig &>(config);
    PipelinePlan plan;
    plan.report_label      = "quality-metrics";
    plan.reporter          = cfg.report.reporter;
    plan.save_run_report   = cfg.report.save_run_report;
    plan.run_report_prefix = "quality-metrics";
    plan.threads           = static_cast<std::size_t>(cfg.runtime.threads);
    plan.success_message   = "Quality-metrics computation completed successfully!";

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
      throw CliUsageError("quality-metrics does not support this source mode");
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
    metrics_task.name = "quality_metrics";
    metrics_task.task = quality_metrics::make_quality_metrics_task(cfg.metrics_config);
    if (needs_parse_task) {
      metrics_task.deps.emplace_back("parse_model");
    }
    metrics_task.thread_safe = true;
    plan.tasks.push_back(std::move(metrics_task));

    PipelineSink sink;
    sink.channel      = std::string(quality_metrics::OutputChannel);
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

const CommandSpec &get_quality_metrics_spec() noexcept {
  static const QualityMetricsSpec spec;
  return spec;
}

} // namespace lahuta::cli
