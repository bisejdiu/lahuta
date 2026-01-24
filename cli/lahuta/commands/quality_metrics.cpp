#include <algorithm>
#include <cctype>
#include <cstdint>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include "analysis/extract/extract_tasks.hpp"
#include "cli/arg_validation.hpp"
#include "cli/extension_utils.hpp"
#include "cli/run_report.hpp"
#include "cli/time_utils.hpp"
#include "commands/quality_metrics.hpp"
#include "commands/reporting.hpp"
#include "logging/logging.hpp"
#include "models/dssp.hpp"
#include "models/plddt.hpp"
#include "pipeline/data_requirements.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/dynamic/sources.hpp"
#include "runtime.hpp"
#include "serialization/json.hpp"
#include "sinks/logging.hpp"
#include "sinks/ndjson.hpp"
#include "utils/span.hpp"

// clang-format off
namespace lahuta::cli {

namespace dyn = pipeline::dynamic;
using pipeline::DataField;
using pipeline::DataFieldSet;

struct QualityMetricsOptions {
  enum class SourceMode { Directory, Vector, FileList, Database };

  SourceMode source_mode = SourceMode::Directory;
  std::string directory_path;
  std::vector<std::string> extensions{".cif", ".cif.gz", ".pdb", ".pdb.gz"};
  bool recursive = false;
  std::vector<std::string> file_vector;
  std::string file_list_path;
  std::string database_path;

  bool is_af2_model = false;
  bool output_stdout = false;
  std::string output_path = "per_protein_metrics.jsonl";
  bool output_override = false;
  bool save_run_report = false;
  bool compute_overlaps = true;
  std::size_t segment_min = 30;
  const PipelineReporter* reporter = nullptr;
  int threads = 8;
  size_t batch_size = 512;
  size_t writer_threads = 1;
};

namespace {

struct GroupSpec {
  std::string name;
  std::uint32_t mask = 0;
  bool segment = false;
  std::string fraction_key;
  std::string segment_count_key;
  std::string max_segment_key;
  std::string long_segment_fraction_key;
};

struct OverlapSpec {
  std::string key;
  std::size_t plddt_index = 0;
  std::size_t dssp_index = 0;
};

struct QualityMetricsConfig {
  std::vector<GroupSpec> plddt_groups;
  std::vector<GroupSpec> dssp_groups;
  std::vector<OverlapSpec> overlaps;
  bool compute_overlaps = true;
  std::size_t segment_min = 30;
};

struct SegmentStats {
  std::size_t segment_count = 0;
  std::size_t max_segment_len = 0;
  std::size_t long_segment_residues = 0;
};

constexpr std::string_view OUTPUT_CHANNEL = "per_protein_metrics";

void initialize_runtime(int num_threads) {
  LahutaRuntime::ensure_initialized(static_cast<std::size_t>(num_threads));
}

std::string normalize_group_name(std::string_view raw) {
  raw = detail::trim_view(raw);

  std::string out;
  out.reserve(raw.size());

  for (unsigned char uc : raw) {
    if (std::isspace(uc)) return {};

    const bool allowed = std::isalnum(uc) || uc == '-' || uc == '_';
    if (!allowed) return {};

    out.push_back(static_cast<char>(std::tolower(uc)));
  }

  return out;
}

std::string normalize_enum_token(std::string_view raw) {
  raw = detail::trim_view(raw);
  std::string out;
  out.reserve(raw.size());
  for (unsigned char uc : raw) {
    if (!std::isalnum(uc)) continue;
    out.push_back(static_cast<char>(std::toupper(uc)));
  }
  return out;
}

bool parse_plddt_token(std::string_view raw, pLDDTCategory& out) {
  const std::string token = normalize_enum_token(raw);

  static constexpr std::pair<std::string_view, pLDDTCategory> pLDDTmap[] = {
    {"VERYHIGH", pLDDTCategory::VeryHigh},
    {"HIGH",     pLDDTCategory::High},
    {"LOW",      pLDDTCategory::Low},
    {"VERYLOW",  pLDDTCategory::VeryLow},
  };

  for (auto [k, v] : pLDDTmap) {
    if (token == k) {
      out = v; return true;
    }
  }
  return false;
}

bool parse_dssp_token(std::string_view raw, DSSPAssignment& out) {
  const std::string token = normalize_enum_token(raw);

  static constexpr std::pair<std::string_view, DSSPAssignment> kMap[] = {
    {"COIL",             DSSPAssignment::Coil},
    {"ALPHAHELIX",       DSSPAssignment::AlphaHelix},
    {"HELIX310",         DSSPAssignment::Helix3_10},
    {"HELIXPI",          DSSPAssignment::HelixPi},
    {"POLYPROLINEHELIX", DSSPAssignment::PolyProlineHelix},
    {"STRAND",           DSSPAssignment::Strand},
    {"TURN",             DSSPAssignment::Turn},
    {"BEND",             DSSPAssignment::Bend},
  };

  for (auto [k, v] : kMap) {
    if (token == k) {
      out = v; return true;
    }
  }
  return false;
}

std::uint32_t plddt_mask_for(pLDDTCategory cat) {
  return 1u << static_cast<std::uint32_t>(cat);
}

std::uint32_t dssp_mask_for(DSSPAssignment cat) {
  return 1u << static_cast<std::uint32_t>(cat);
}

bool add_group_definition(std::string_view raw,
                          bool is_plddt,
                          std::unordered_set<std::string>& seen_names,
                          std::vector<GroupSpec>& groups) {
  raw = detail::trim_view(raw);
  if (raw.empty()) {
    Logger::get_logger()->error("Group definition cannot be empty.");
    return false;
  }

  const auto eq_pos = raw.find('=');
  if (eq_pos == std::string_view::npos) {
    Logger::get_logger()->error("Group definition '{}' is missing '=' (expected name=values).", raw);
    return false;
  }

  const std::string name = normalize_group_name(raw.substr(0, eq_pos));
  if (name.empty()) {
    Logger::get_logger()->error("Group name '{}' is invalid (use letters, digits, '-' or '_').", raw.substr(0, eq_pos));
    return false;
  }
  if (!seen_names.insert(name).second) {
    Logger::get_logger()->error("Duplicate group name '{}'.", name);
    return false;
  }

  std::vector<std::string> tokens;
  detail::split_argument_list(raw.substr(eq_pos + 1), false, tokens);
  if (tokens.empty()) {
    Logger::get_logger()->error("Group '{}' has no values.", name);
    return false;
  }

  std::uint32_t mask = 0;
  for (const auto& token : tokens) {
    if (is_plddt) {
      pLDDTCategory cat{};
      if (!parse_plddt_token(token, cat)) {
        Logger::get_logger()->error("Invalid pLDDT token '{}' in group '{}'.", token, name);
        return false;
      }
      mask |= plddt_mask_for(cat);
    } else {
      DSSPAssignment cat{};
      if (!parse_dssp_token(token, cat)) {
        Logger::get_logger()->error("Invalid DSSP token '{}' in group '{}'.", token, name);
        return false;
      }
      mask |= dssp_mask_for(cat);
    }
  }

  if (mask == 0) {
    Logger::get_logger()->error("Group '{}' resulted in an empty mask.", name);
    return false;
  }

  GroupSpec group;
  group.name = name;
  group.mask = mask;
  groups.push_back(std::move(group));
  return true;
}

void add_default_plddt_groups(std::vector<GroupSpec>& groups) {
  GroupSpec poor;
  poor.name = "poor";
  poor.mask = plddt_mask_for(pLDDTCategory::VeryLow);
  groups.push_back(std::move(poor));

  GroupSpec low_or_poor;
  low_or_poor.name = "low_or_poor";
  low_or_poor.mask = plddt_mask_for(pLDDTCategory::Low) | plddt_mask_for(pLDDTCategory::VeryLow);
  groups.push_back(std::move(low_or_poor));
}

void add_default_dssp_groups(std::vector<GroupSpec>& groups) {
  GroupSpec coil;
  coil.name = "coil";
  coil.mask = dssp_mask_for(DSSPAssignment::Coil)
              | dssp_mask_for(DSSPAssignment::Turn)
              | dssp_mask_for(DSSPAssignment::Bend);
  groups.push_back(std::move(coil));

  GroupSpec helix;
  helix.name = "helix";
  helix.mask = dssp_mask_for(DSSPAssignment::AlphaHelix)
               | dssp_mask_for(DSSPAssignment::Helix3_10)
               | dssp_mask_for(DSSPAssignment::HelixPi)
               | dssp_mask_for(DSSPAssignment::PolyProlineHelix);
  groups.push_back(std::move(helix));

  GroupSpec strand;
  strand.name = "strand";
  strand.mask = dssp_mask_for(DSSPAssignment::Strand);
  groups.push_back(std::move(strand));
}

void build_plddt_group_keys(GroupSpec& group) {
  group.fraction_key = "plddt_" + group.name + "_fraction";
  if (group.segment) {
    group.segment_count_key = "plddt_" + group.name + "_segment_count";
    group.max_segment_key = "plddt_" + group.name + "_max_segment_len";
    group.long_segment_fraction_key = "plddt_" + group.name + "_long_segment_fraction";
  }
}

void build_dssp_group_keys(GroupSpec& group) {
  group.fraction_key = "dssp_" + group.name + "_fraction";
}

std::string build_overlap_key(std::string_view plddt, std::string_view dssp) {
  std::string key;
  key.reserve(plddt.size() + dssp.size() + 24);
  key.append("overlap_");
  key.append(plddt);
  key.append("__");
  key.append(dssp);
  key.append("_fraction");
  return key;
}

SegmentStats compute_segment_stats(span<const pLDDTCategory> plddt,
                                   std::uint32_t mask,
                                   std::size_t long_threshold) {
  SegmentStats stats;
  std::size_t current = 0;

  for (const auto& cat : plddt) {
    const bool in_group = (mask & plddt_mask_for(cat)) != 0;
    if (in_group) {
      current += 1;
      continue;
    }
    if (current > 0) {
      stats.segment_count += 1;
      stats.max_segment_len = std::max(stats.max_segment_len, current);
      if (current >= long_threshold) {
        stats.long_segment_residues += current;
      }
      current = 0;
    }
  }

  if (current > 0) {
    stats.segment_count += 1;
    stats.max_segment_len = std::max(stats.max_segment_len, current);
    if (current >= long_threshold) {
      stats.long_segment_residues += current;
    }
  }

  return stats;
}

class QualityMetricsTask final : public dyn::ITask {
public:
  explicit QualityMetricsTask(std::shared_ptr<const QualityMetricsConfig> config)
      : config_(std::move(config)) {}

  dyn::TaskResult run(const std::string& item_path, dyn::TaskContext& ctx) override {
    auto payload = ctx.model_payload();
    std::shared_ptr<const ModelParserResult> parsed;

    span<const pLDDTCategory> plddt_span;
    span<const DSSPAssignment> dssp_span;

    if (payload) {
      if (payload->plddts && !payload->plddts->empty()) {
        plddt_span = span(*payload->plddts);
      }

      if (payload->dssp && !payload->dssp->empty()) {
        dssp_span = span(*payload->dssp);
      }
    } else {
      parsed = analysis::extract::get_cached_model_parser_result(ctx);
      if (parsed) {
        if (!parsed->plddt_per_residue.empty()) {
          plddt_span = span(std::as_const(parsed->plddt_per_residue));
        }
        if (!parsed->dssp_per_residue.empty()) {
          dssp_span = span(std::as_const(parsed->dssp_per_residue));
        }
      }
    }

    if (!config_->plddt_groups.empty() && plddt_span.empty()) {
      Logger::get_logger()->warn("[quality-metrics] Missing pLDDT data for '{}'", item_path);
      return {};
    }
    if (!config_->dssp_groups.empty() && dssp_span.empty()) {
      Logger::get_logger()->warn("[quality-metrics] Missing DSSP data for '{}'", item_path);
      return {};
    }

    std::size_t length = 0;
    if      (!plddt_span.empty()) length = plddt_span.size();
    else if (!dssp_span .empty()) length = dssp_span.size();

    if (!plddt_span.empty() && !dssp_span.empty() && plddt_span.size() != dssp_span.size()) {
      Logger::get_logger()->warn("[quality-metrics] Length mismatch for '{}': pLDDT={} DSSP={}",
                                 item_path, plddt_span.size(), dssp_span.size());
      return {};
    }

    if (length == 0) {
      Logger::get_logger()->warn("[quality-metrics] Empty model payload for '{}'", item_path);
      return {};
    }

    std::vector<std::size_t> plddt_counts(config_->plddt_groups.size(), 0);
    std::vector<std::size_t> dssp_counts (config_-> dssp_groups.size(), 0);
    std::vector<std::size_t> overlap_counts;
    if (config_->compute_overlaps && !config_->plddt_groups.empty() && !config_->dssp_groups.empty()) {
      overlap_counts.assign(config_->plddt_groups.size() * config_->dssp_groups.size(), 0);
    }

    std::vector<std::size_t> plddt_matches;
    std::vector<std::size_t> dssp_matches;
    plddt_matches.reserve(config_->plddt_groups.size());
    dssp_matches.reserve(config_->dssp_groups.size());

    for (std::size_t idx = 0; idx < length; ++idx) {
      plddt_matches.clear();
      dssp_matches.clear();

      if (!plddt_span.empty()) {
        const auto cat = plddt_span[idx];
        const std::uint32_t bit = plddt_mask_for(cat);
        for (std::size_t i = 0; i < config_->plddt_groups.size(); ++i) {
          if (config_->plddt_groups[i].mask & bit) {
            plddt_counts[i] += 1;
            if (!overlap_counts.empty()) plddt_matches.push_back(i);
          }
        }
      }

      if (!dssp_span.empty()) {
        const auto cat = dssp_span[idx];
        const std::uint32_t bit = dssp_mask_for(cat);
        for (std::size_t i = 0; i < config_->dssp_groups.size(); ++i) {
          if (config_->dssp_groups[i].mask & bit) {
            dssp_counts[i] += 1;
            if (!overlap_counts.empty()) dssp_matches.push_back(i);
          }
        }
      }

      if (!overlap_counts.empty()) {
        for (const auto p_idx : plddt_matches) {
          const std::size_t base = p_idx * config_->dssp_groups.size();
          for (const auto d_idx : dssp_matches) {
            overlap_counts[base + d_idx] += 1;
          }
        }
      }
    }

    std::vector<SegmentStats> segment_stats(config_->plddt_groups.size());
    if (!plddt_span.empty()) {
      for (std::size_t i = 0; i < config_->plddt_groups.size(); ++i) {
        if (!config_->plddt_groups[i].segment) continue;
        segment_stats[i] = compute_segment_stats(plddt_span, config_->plddt_groups[i].mask, config_->segment_min);
      }
    }

    std::string organism = "Unknown";
    std::string taxonomy_id = "N/A";
    if (payload && payload->metadata) {
      organism = payload->metadata->organism_scientific.empty() ? "Unknown" : payload->metadata->organism_scientific;
      taxonomy_id = payload->metadata->ncbi_taxonomy_id.empty() ? "N/A" : payload->metadata->ncbi_taxonomy_id;
    } else if (parsed) {
      organism = parsed->metadata.organism_scientific.empty() ? "Unknown" : parsed->metadata.organism_scientific;
      taxonomy_id = parsed->metadata.ncbi_taxonomy_id.empty() ? "N/A" : parsed->metadata.ncbi_taxonomy_id;
    }

    JsonBuilder json(512);
    json.key("model").value(item_path)
        .key("length").value(static_cast<std::uint64_t>(length))
        .key("organism").value(organism)
        .key("taxonomy_id").value(taxonomy_id);

    for (std::size_t i = 0; i < config_->plddt_groups.size(); ++i) {
      const auto& group = config_->plddt_groups[i];
      const double fraction = static_cast<double>(plddt_counts[i]) / static_cast<double>(length);
      json.key(group.fraction_key).value(fraction);

      if (group.segment) {
        const auto& stats = segment_stats[i];
        const double long_fraction = static_cast<double>(stats.long_segment_residues) / static_cast<double>(length);
        json.key(group.segment_count_key).value(static_cast<std::uint64_t>(stats.segment_count))
            .key(group.max_segment_key).value(static_cast<std::uint64_t>(stats.max_segment_len))
            .key(group.long_segment_fraction_key).value(long_fraction);
      }
    }

    for (std::size_t i = 0; i < config_->dssp_groups.size(); ++i) {
      const auto& group = config_->dssp_groups[i];
      const double fraction = static_cast<double>(dssp_counts[i]) / static_cast<double>(length);
      json.key(group.fraction_key).value(fraction);
    }

    if (!overlap_counts.empty()) {
      for (const auto& overlap : config_->overlaps) {
        const std::size_t index = overlap.plddt_index * config_->dssp_groups.size() + overlap.dssp_index;
        const double fraction = static_cast<double>(overlap_counts[index]) / static_cast<double>(length);
        json.key(overlap.key).value(fraction);
      }
    }

    dyn::TaskResult result;
    result.ok = true;
    result.emits.push_back(dyn::Emission{std::string(OUTPUT_CHANNEL), json.str()});
    return result;
  }

  DataFieldSet data_requirements() const override {
    DataFieldSet fields = DataFieldSet::none();
    if (!config_->plddt_groups.empty()) fields |= DataField::Plddt;
    if (!config_->dssp_groups.empty()) fields |= DataField::Dssp;
    fields |= DataField::Metadata;
    return fields;
  }

private:
  std::shared_ptr<const QualityMetricsConfig> config_;
};

} // namespace

namespace quality_metrics_opts {
const option::Descriptor usage[] = {
  {QualityMetricsOptionIndex::Unknown, 0, "", "", validate::Unknown,
   "Usage: lahuta quality-metrics [options]\n\n"
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
   "    --output metrics.jsonl\n\n"
   "Input Options (choose one):"},
  {QualityMetricsOptionIndex::Help, 0, "h", "help", option::Arg::None,
   "  --help, -h                   \tPrint this help message and exit."},
  {QualityMetricsOptionIndex::PlddtGroup, 0, "", "plddt-group", validate::Required,
   "  --plddt-group <name=values>  \tAdd a pLDDT group (repeatable)."},
  {QualityMetricsOptionIndex::DsspGroup, 0, "", "dssp-group", validate::Required,
   "  --dssp-group <name=values>   \tAdd a DSSP group (repeatable)."},
  {QualityMetricsOptionIndex::SegmentGroup, 0, "", "segment-group", validate::Required,
   "  --segment-group <name>       \tEnable segment metrics for a pLDDT group (repeatable)."},
  {QualityMetricsOptionIndex::SegmentMin, 0, "", "segment-min", validate::Required,
   "  --segment-min <N>            \tMinimum residues for long-segment fraction (default: 30)."},
  {QualityMetricsOptionIndex::NoOverlap, 0, "", "no-overlap", option::Arg::None,
   "  --no-overlap                 \tDisable pLDDT x DSSP overlap metrics."},
  {QualityMetricsOptionIndex::SourceDatabase, 0, "", "database", validate::Required,
   "  --database <path>            \tProcess structures from database."},
  {QualityMetricsOptionIndex::SourceDirectory, 0, "d", "directory", validate::Required,
   "  --directory, -d <path>       \tProcess all files in directory."},
  {QualityMetricsOptionIndex::SourceVector, 0, "f", "files", validate::Required,
   "  --files, -f <file1,file2>    \tProcess specific files (comma-separated or repeat -f)."},
  {QualityMetricsOptionIndex::SourceFileList, 0, "l", "file-list", validate::Required,
   "  --file-list, -l <path>       \tProcess files listed in text file (one per line)."},
  {0, 0, "", "", option::Arg::None,
   "\nDirectory Options:"},
  {QualityMetricsOptionIndex::Extension, 0, "e", "extension", validate::Required,
   "  --extension, -e <ext>        \tFile extension(s) for directory mode. Repeat or comma-separate values (default: .cif, .cif.gz, .pdb, .pdb.gz)."},
  {QualityMetricsOptionIndex::Recursive, 0, "r", "recursive", option::Arg::None,
   "  --recursive, -r              \tRecursively search subdirectories."},
  {0, 0, "", "", option::Arg::None,
   "\nModel Options:"},
  {QualityMetricsOptionIndex::IsAf2Model, 0, "", "is_af2_model", option::Arg::None,
   "  --is_af2_model               \tRequired for file-based sources; inputs are AlphaFold2 models (AF2-like mmCIF)."},
  {0, 0, "", "", option::Arg::None,
   "\nOutput Options:"},
  {QualityMetricsOptionIndex::Output, 0, "o", "output", validate::Required,
   "  --output, -o <path>          \tWrite NDJSON to file (default: per_protein_metrics.jsonl). Use '-' for stdout."},
  {QualityMetricsOptionIndex::Reporter, 0, "", "reporter", validate::Required,
   "  --reporter <name>            \tSelect pipeline reporter. See help for names."},
  {QualityMetricsOptionIndex::SaveRunReport, 0, "", "save-run-report", option::Arg::None,
   "  --save-run-report            \tSave run statistics to a JSON file."},
  {0, 0, "", "", option::Arg::None,
   "\nRuntime Options:"},
  {QualityMetricsOptionIndex::Threads, 0, "t", "threads", validate::Required,
   "  --threads, -t <num>          \tNumber of threads to use (default: 8)."},
  {QualityMetricsOptionIndex::BatchSize, 0, "b", "batch-size", validate::Required,
   "  --batch-size, -b <size>      \tBatch size for processing (default: 512)."},
  {QualityMetricsOptionIndex::WriterThreads, 0, "", "writer-threads", validate::Required,
   "  --writer-threads <num>       \tNumber of writer threads per sink (default: 1)."},
  {0, 0, 0, 0, 0, 0}
};
} // namespace quality_metrics_opts

[[nodiscard]] std::unique_ptr<CliCommand> QualityMetricsCommand::create() {
  return std::unique_ptr<CliCommand>(new QualityMetricsCommand());
}

int QualityMetricsCommand::run(int argc, char* argv[]) {
  option::Stats stats(true, quality_metrics_opts::usage, argc, const_cast<const char**>(argv));
  std::vector<option::Option> options(stats.options_max);
  std::vector<option::Option> buffer (stats.buffer_max);
  option::Parser parse(true, quality_metrics_opts::usage, argc, const_cast<const char**>(argv), options.data(), buffer.data());

  if (parse.error()) return 1;

  if (options[quality_metrics_opts::QualityMetricsOptionIndex::Help]) {
    option::printUsage(std::cout, quality_metrics_opts::usage);
    return 0;
  }

  try {
    if (parse.nonOptionsCount() > 0) {
      Logger::get_logger()->error("Unexpected positional argument '{}'.", parse.nonOption(0));
      option::printUsage(std::cerr, quality_metrics_opts::usage);
      return 1;
    }

    QualityMetricsOptions cli;
    const auto default_sink_cfg = dyn::get_default_backpressure_config();
    cli.writer_threads = default_sink_cfg.writer_threads;
    cli.reporter = &default_pipeline_reporter();

    std::unordered_set<std::string> plddt_group_names;
    std::unordered_set<std::string> dssp_group_names;
    std::vector<GroupSpec> plddt_groups;
    std::vector<GroupSpec> dssp_groups;

    for (const option::Option* opt = &options[quality_metrics_opts::QualityMetricsOptionIndex::PlddtGroup];
         opt != nullptr;
         opt = opt->next()) {
      if (!opt->arg) continue;
      if (!add_group_definition(opt->arg, true, plddt_group_names, plddt_groups)) {
        return 1;
      }
    }

    for (const option::Option* opt = &options[quality_metrics_opts::QualityMetricsOptionIndex::DsspGroup];
         opt != nullptr;
         opt = opt->next()) {
      if (!opt->arg) continue;
      if (!add_group_definition(opt->arg, false, dssp_group_names, dssp_groups)) {
        return 1;
      }
    }

    bool used_default_plddt = false;
    if (plddt_groups.empty()) {
      add_default_plddt_groups(plddt_groups);
      used_default_plddt = true;
    }
    if (dssp_groups.empty()) {
      add_default_dssp_groups(dssp_groups);
    }

    std::unordered_set<std::string> segment_group_names;
    for (const option::Option* opt = &options[quality_metrics_opts::QualityMetricsOptionIndex::SegmentGroup];
         opt != nullptr;
         opt = opt->next()) {
      if (!opt->arg) continue;
      std::vector<std::string> tokens;
      detail::split_argument_list(opt->arg, false, tokens);
      for (const auto& token : tokens) {
        const std::string name = normalize_group_name(token);
        if (name.empty()) {
          Logger::get_logger()->error("Invalid segment group name '{}'.", token);
          return 1;
        }
        segment_group_names.insert(name);
      }
    }

    if (segment_group_names.empty() && used_default_plddt) {
      segment_group_names.insert("poor");
    }

    if (!segment_group_names.empty()) {
      std::unordered_set<std::string> available;
      for (const auto& group : plddt_groups) {
        available.insert(group.name);
      }
      for (const auto& name : segment_group_names) {
        if (available.count(name) == 0) {
          Logger::get_logger()->error("Segment group '{}' does not match any pLDDT group.", name);
          return 1;
        }
      }
      for (auto& group : plddt_groups) {
        group.segment = segment_group_names.count(group.name) > 0;
      }
    }

    for (auto& group : plddt_groups) build_plddt_group_keys(group);
    for (auto& group : dssp_groups)  build_dssp_group_keys(group);

    cli.compute_overlaps = options[quality_metrics_opts::QualityMetricsOptionIndex::NoOverlap] ? false : true;

    if (options[quality_metrics_opts::QualityMetricsOptionIndex::SegmentMin]) {
      cli.segment_min = std::stoull(options[quality_metrics_opts::QualityMetricsOptionIndex::SegmentMin].arg);
      if (cli.segment_min == 0) {
        Logger::get_logger()->error("--segment-min must be positive");
        return 1;
      }
    }

    if (options[quality_metrics_opts::QualityMetricsOptionIndex::Output]) {
      std::string_view output_arg = options[quality_metrics_opts::QualityMetricsOptionIndex::Output].arg
                                      ? options[quality_metrics_opts::QualityMetricsOptionIndex::Output].arg
                                      : std::string_view{};
      if (output_arg.empty()) {
        Logger::get_logger()->error("--output requires a value.");
        return 1;
      }
      if (output_arg == "-") {
        cli.output_stdout = true;
        cli.output_path = "-";
      } else {
        cli.output_path = std::string(output_arg);
        cli.output_override = true;
      }
    }

    if (options[quality_metrics_opts::QualityMetricsOptionIndex::Reporter]) {
      std::string_view name = options[quality_metrics_opts::QualityMetricsOptionIndex::Reporter].arg
                                ? options[quality_metrics_opts::QualityMetricsOptionIndex::Reporter].arg
                                : std::string_view{};
      if (name.empty()) {
        Logger::get_logger()->error("--reporter requires a value.");
        return 1;
      }
      if (const auto* rep = find_pipeline_reporter(name)) {
        cli.reporter = rep;
      } else {
        Logger::get_logger()->error("Unknown reporter '{}' ", name);
        return 1;
      }
    }

    cli.save_run_report = options[quality_metrics_opts::QualityMetricsOptionIndex::SaveRunReport] ? true : false;

    // Parse source options
    int source_count = 0;
    if (options[quality_metrics_opts::QualityMetricsOptionIndex::SourceDatabase]) {
      cli.source_mode = QualityMetricsOptions::SourceMode::Database;
      cli.database_path = options[quality_metrics_opts::QualityMetricsOptionIndex::SourceDatabase].arg;
      source_count++;
    }
    if (options[quality_metrics_opts::QualityMetricsOptionIndex::SourceDirectory]) {
      cli.source_mode = QualityMetricsOptions::SourceMode::Directory;
      cli.directory_path = options[quality_metrics_opts::QualityMetricsOptionIndex::SourceDirectory].arg;
      source_count++;
    }
    if (options[quality_metrics_opts::QualityMetricsOptionIndex::SourceVector]) {
      cli.source_mode = QualityMetricsOptions::SourceMode::Vector;
      cli.file_vector.clear();
      for (const option::Option* opt = &options[quality_metrics_opts::QualityMetricsOptionIndex::SourceVector];
           opt != nullptr;
           opt = opt->next()) {
        if (opt->arg) parse_file_argument(opt->arg, cli.file_vector);
      }
      source_count++;
    }
    if (options[quality_metrics_opts::QualityMetricsOptionIndex::SourceFileList]) {
      cli.source_mode = QualityMetricsOptions::SourceMode::FileList;
      cli.file_list_path = options[quality_metrics_opts::QualityMetricsOptionIndex::SourceFileList].arg;
      source_count++;
    }

    if (source_count == 0) {
      Logger::get_logger()->error("Must specify exactly one source option: --database, --directory, --files, or --file-list");
      return 1;
    }
    if (source_count > 1) {
      Logger::get_logger()->error("Cannot specify multiple source options");
      return 1;
    }

    if (options[quality_metrics_opts::QualityMetricsOptionIndex::Extension]) {
      cli.extensions.clear();
      for (const option::Option* opt = &options[quality_metrics_opts::QualityMetricsOptionIndex::Extension]; opt != nullptr; opt = opt->next()) {
        parse_extension_argument(opt->arg ? opt->arg : "", cli.extensions);
      }
    }
    if (cli.extensions.empty()) cli.extensions.emplace_back();

    cli.recursive    = options[quality_metrics_opts::QualityMetricsOptionIndex::Recursive]  ? true : false;
    cli.is_af2_model = options[quality_metrics_opts::QualityMetricsOptionIndex::IsAf2Model] ? true : false;

    if (cli.source_mode == QualityMetricsOptions::SourceMode::Database) {
      cli.is_af2_model = true;
    }

    if (cli.source_mode != QualityMetricsOptions::SourceMode::Database && !cli.is_af2_model) {
      Logger::get_logger()->error(
          "quality-metrics expects AlphaFold2 model inputs. For file-based sources, pass --is_af2_model (or use --database).");
      option::printUsage(std::cerr, quality_metrics_opts::usage);
      return 1;
    }

    if (options[quality_metrics_opts::QualityMetricsOptionIndex::Threads]) {
      cli.threads = std::stoi(options[quality_metrics_opts::QualityMetricsOptionIndex::Threads].arg);
      if (cli.threads <= 0) {
        Logger::get_logger()->error("Threads must be positive");
        return 1;
      }
    }

    if (options[quality_metrics_opts::QualityMetricsOptionIndex::BatchSize]) {
      cli.batch_size = std::stoull(options[quality_metrics_opts::QualityMetricsOptionIndex::BatchSize].arg);
      if (cli.batch_size == 0) {
        Logger::get_logger()->error("Batch size must be positive");
        return 1;
      }
    }

    if (options[quality_metrics_opts::QualityMetricsOptionIndex::WriterThreads]) {
      cli.writer_threads = std::stoull(options[quality_metrics_opts::QualityMetricsOptionIndex::WriterThreads].arg);
      if (cli.writer_threads == 0) {
        Logger::get_logger()->error("Writer threads must be positive");
        return 1;
      }
    }

    initialize_runtime(cli.threads);

    auto source = std::unique_ptr<sources::IDescriptor>{};
    switch (cli.source_mode) {
      case QualityMetricsOptions::SourceMode::Database:
        source = dyn::sources_factory::from_lmdb(
            cli.database_path,
            std::string{},
            cli.batch_size,
            {static_cast<std::size_t>(cli.threads) + 1});
        break;
      case QualityMetricsOptions::SourceMode::Directory:
        source = dyn::sources_factory::from_directory(
            cli.directory_path,
            cli.extensions,
            cli.recursive,
            cli.batch_size);
        break;
      case QualityMetricsOptions::SourceMode::Vector:
        source = dyn::sources_factory::from_vector(cli.file_vector);
        break;
      case QualityMetricsOptions::SourceMode::FileList:
        source = dyn::sources_factory::from_filelist(cli.file_list_path);
        break;
    }

    dyn::StageManager mgr(std::move(source));
    auto sink_cfg = dyn::get_default_backpressure_config();
    sink_cfg.writer_threads = cli.writer_threads;

    const bool needs_parse_task = cli.source_mode != QualityMetricsOptions::SourceMode::Database;
    if (needs_parse_task) {
      auto parse_task = std::make_shared<analysis::extract::ModelParseTask>();
      mgr.add_task("parse_model", {}, std::move(parse_task), /*thread_safe=*/true);
    }

    auto config = std::make_shared<QualityMetricsConfig>();
    config->plddt_groups = std::move(plddt_groups);
    config->dssp_groups = std::move(dssp_groups);
    config->compute_overlaps = cli.compute_overlaps;
    config->segment_min = cli.segment_min;
    if (config->compute_overlaps && !config->plddt_groups.empty() && !config->dssp_groups.empty()) {
      config->overlaps.reserve(config->plddt_groups.size() * config->dssp_groups.size());
      for (std::size_t i = 0; i < config->plddt_groups.size(); ++i) {
        for (std::size_t j = 0; j < config->dssp_groups.size(); ++j) {
          OverlapSpec overlap;
          overlap.plddt_index = i;
          overlap.dssp_index = j;
          overlap.key = build_overlap_key(config->plddt_groups[i].name, config->dssp_groups[j].name);
          config->overlaps.push_back(std::move(overlap));
        }
      }
    }

    auto task = std::make_shared<QualityMetricsTask>(std::move(config));
    std::vector<std::string> deps;
    if (needs_parse_task) deps.emplace_back("parse_model");
    mgr.add_task("quality_metrics", std::move(deps), std::move(task), /*thread_safe=*/true);

    std::shared_ptr<dyn::LoggingSink> stdout_sink;
    if (cli.output_stdout) {
      stdout_sink = std::make_shared<dyn::LoggingSink>();
      Logger::get_logger()->info("Quality-metrics output sink -> stdout (-)");
      mgr.connect_sink(std::string(OUTPUT_CHANNEL), stdout_sink, sink_cfg);
    } else {
      Logger::get_logger()->info("Quality-metrics output sink -> file: {}", cli.output_path);
      mgr.connect_sink(std::string(OUTPUT_CHANNEL),
                       std::make_shared<dyn::NdjsonFileSink>(cli.output_path),
                       sink_cfg);
    }

    mgr.compile();
    auto progress = attach_progress_observer(mgr);
    const auto report = mgr.run(static_cast<std::size_t>(cli.threads));
    if (progress) progress->finish();
    const auto* reporter = cli.reporter ? cli.reporter : &default_pipeline_reporter();
    reporter->emit("quality-metrics", report);
    if (cli.save_run_report) {
      const std::string report_path = make_report_path("quality-metrics", report.run_token, current_timestamp_string());
      if (!write_run_report_json(report_path, report)) {
        throw std::runtime_error("Failed to persist RunReport JSON");
      }
      Logger::get_logger()->info("Run report saved to {}", report_path);
    }

    return 0;
  } catch (const std::exception& e) {
    Logger::get_logger()->error("Error: {}", e.what());
    return 1;
  }
}

} // namespace lahuta::cli
