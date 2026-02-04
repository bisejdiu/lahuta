/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string dst;
 *   std::for_each(parts.begin(), parts.end(), [&dst](std::string_view p) { dst += p; });
 *   return dst;
 * }();
 *
 */

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "analysis/system/model_parse_task.hpp"
#include "logging/logging.hpp"
#include "models/dssp.hpp"
#include "models/plddt.hpp"
#include "parsing/extension_utils.hpp"
#include "pipeline/data/data_requirements.hpp"
#include "serialization/json.hpp"
#include "tasks/quality_metrics_task.hpp"
#include "utils/span.hpp"

namespace lahuta::cli::quality_metrics {
namespace A = lahuta::analysis;
namespace P = lahuta::pipeline;
namespace {

struct SegmentStats {
  std::size_t segment_count         = 0;
  std::size_t max_segment_len       = 0;
  std::size_t long_segment_residues = 0;
};

std::string normalize_enum_token(std::string_view raw) {
  raw = lahuta::cli::detail::trim_view(raw);
  std::string out;
  out.reserve(raw.size());
  for (unsigned char uc : raw) {
    if (!std::isalnum(uc)) continue;
    out.push_back(static_cast<char>(std::toupper(uc)));
  }
  return out;
}

bool parse_plddt_token(std::string_view raw, pLDDTCategory &out) {
  const std::string token = normalize_enum_token(raw);

  static constexpr std::pair<std::string_view, pLDDTCategory> plddt_map[] = {
      {"VERYHIGH", pLDDTCategory::VeryHigh},
      {"HIGH",     pLDDTCategory::High    },
      {"LOW",      pLDDTCategory::Low     },
      {"VERYLOW",  pLDDTCategory::VeryLow },
  };

  for (auto [k, v] : plddt_map) {
    if (token == k) {
      out = v;
      return true;
    }
  }
  return false;
}

bool parse_dssp_token(std::string_view raw, DSSPAssignment &out) {
  const std::string token = normalize_enum_token(raw);

  static constexpr std::pair<std::string_view, DSSPAssignment> Map[] = {
      {"COIL",             DSSPAssignment::Coil            },
      {"ALPHAHELIX",       DSSPAssignment::AlphaHelix      },
      {"HELIX310",         DSSPAssignment::Helix3_10       },
      {"HELIXPI",          DSSPAssignment::HelixPi         },
      {"POLYPROLINEHELIX", DSSPAssignment::PolyProlineHelix},
      {"STRAND",           DSSPAssignment::Strand          },
      {"TURN",             DSSPAssignment::Turn            },
      {"BEND",             DSSPAssignment::Bend            },
  };

  for (auto [k, v] : Map) {
    if (token == k) {
      out = v;
      return true;
    }
  }
  return false;
}

std::uint32_t plddt_mask_for(pLDDTCategory cat) { return 1u << static_cast<std::uint32_t>(cat); }

std::uint32_t dssp_mask_for(DSSPAssignment cat) { return 1u << static_cast<std::uint32_t>(cat); }

void build_plddt_group_keys(GroupSpec &group) {
  group.fraction_key = "plddt_" + group.name + "_fraction";
  if (group.segment) {
    group.segment_count_key         = "plddt_" + group.name + "_segment_count";
    group.max_segment_key           = "plddt_" + group.name + "_max_segment_len";
    group.long_segment_fraction_key = "plddt_" + group.name + "_long_segment_fraction";
  }
}

void build_dssp_group_keys(GroupSpec &group) { group.fraction_key = "dssp_" + group.name + "_fraction"; }

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

SegmentStats compute_segment_stats(span<const pLDDTCategory> plddt, std::uint32_t mask,
                                   std::size_t long_threshold) {
  SegmentStats stats;
  std::size_t current = 0;

  for (const auto &cat : plddt) {
    const bool in_group = (mask & plddt_mask_for(cat)) != 0;
    if (in_group) {
      current += 1;
      continue;
    }
    if (current > 0) {
      stats.segment_count   += 1;
      stats.max_segment_len  = std::max(stats.max_segment_len, current);
      if (current >= long_threshold) {
        stats.long_segment_residues += current;
      }
      current = 0;
    }
  }

  if (current > 0) {
    stats.segment_count   += 1;
    stats.max_segment_len  = std::max(stats.max_segment_len, current);
    if (current >= long_threshold) {
      stats.long_segment_residues += current;
    }
  }

  return stats;
}

class QualityMetricsTask final : public P::ITask {
public:
  explicit QualityMetricsTask(std::shared_ptr<const QualityMetricsConfig> config)
      : config_(std::move(config)) {}

  P::TaskResult run(const std::string &item_path, P::TaskContext &ctx) override {
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
      parsed = A::get_parsed_model_result(ctx);
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
      Logger::get_logger()->warn("[quality-metrics:input] Missing pLDDT data for '{}'", item_path);
      return {};
    }
    if (!config_->dssp_groups.empty() && dssp_span.empty()) {
      Logger::get_logger()->warn("[quality-metrics:input] Missing DSSP data for '{}'", item_path);
      return {};
    }

    std::size_t length = 0;
    if (!plddt_span.empty())
      length = plddt_span.size();
    else if (!dssp_span.empty())
      length = dssp_span.size();

    if (!plddt_span.empty() && !dssp_span.empty() && plddt_span.size() != dssp_span.size()) {
      Logger::get_logger()->warn("[quality-metrics:input] Length mismatch for '{}': pLDDT={} DSSP={}",
                                 item_path,
                                 plddt_span.size(),
                                 dssp_span.size());
      return {};
    }

    if (length == 0) {
      Logger::get_logger()->warn("[quality-metrics:input] Empty model payload for '{}'", item_path);
      return {};
    }

    std::vector<std::size_t> plddt_counts(config_->plddt_groups.size(), 0);
    std::vector<std::size_t> dssp_counts(config_->dssp_groups.size(), 0);
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
        const auto cat          = plddt_span[idx];
        const std::uint32_t bit = plddt_mask_for(cat);
        for (std::size_t i = 0; i < config_->plddt_groups.size(); ++i) {
          if (config_->plddt_groups[i].mask & bit) {
            plddt_counts[i] += 1;
            if (!overlap_counts.empty()) plddt_matches.push_back(i);
          }
        }
      }

      if (!dssp_span.empty()) {
        const auto cat          = dssp_span[idx];
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
        segment_stats[i] = compute_segment_stats(plddt_span,
                                                 config_->plddt_groups[i].mask,
                                                 config_->segment_min);
      }
    }

    std::string organism    = "Unknown";
    std::string taxonomy_id = "N/A";
    if (payload && payload->metadata) {
      organism    = payload->metadata->organism_scientific.empty() ? "Unknown"
                                                                   : payload->metadata->organism_scientific;
      taxonomy_id = payload->metadata->ncbi_taxonomy_id.empty() ? "N/A" : payload->metadata->ncbi_taxonomy_id;
    } else if (parsed) {
      organism    = parsed->metadata.organism_scientific.empty() ? "Unknown"
                                                                 : parsed->metadata.organism_scientific;
      taxonomy_id = parsed->metadata.ncbi_taxonomy_id.empty() ? "N/A" : parsed->metadata.ncbi_taxonomy_id;
    }

    JsonBuilder json(512);
    json.key("model")
        .value(item_path)
        .key("length")
        .value(static_cast<std::uint64_t>(length))
        .key("organism")
        .value(organism)
        .key("taxonomy_id")
        .value(taxonomy_id);

    for (std::size_t i = 0; i < config_->plddt_groups.size(); ++i) {
      const auto &group     = config_->plddt_groups[i];
      const double fraction = static_cast<double>(plddt_counts[i]) / static_cast<double>(length);
      json.key(group.fraction_key).value(fraction);

      if (group.segment) {
        const auto &stats          = segment_stats[i];
        const double long_fraction = static_cast<double>(stats.long_segment_residues) /
                                     static_cast<double>(length);
        json.key(group.segment_count_key)
            .value(static_cast<std::uint64_t>(stats.segment_count))
            .key(group.max_segment_key)
            .value(static_cast<std::uint64_t>(stats.max_segment_len))
            .key(group.long_segment_fraction_key)
            .value(long_fraction);
      }
    }

    for (std::size_t i = 0; i < config_->dssp_groups.size(); ++i) {
      const auto &group     = config_->dssp_groups[i];
      const double fraction = static_cast<double>(dssp_counts[i]) / static_cast<double>(length);
      json.key(group.fraction_key).value(fraction);
    }

    if (!overlap_counts.empty()) {
      for (const auto &overlap : config_->overlaps) {
        const std::size_t index = overlap.plddt_index * config_->dssp_groups.size() + overlap.dssp_index;
        const double fraction   = static_cast<double>(overlap_counts[index]) / static_cast<double>(length);
        json.key(overlap.key).value(fraction);
      }
    }

    P::TaskResult result;
    result.ok = true;
    result.emits.push_back(P::Emission{std::string(OutputChannel), json.str()});
    return result;
  }

  P::DataFieldSet data_requirements() const override {
    P::DataFieldSet fields = P::DataFieldSet::none();
    if (!config_->plddt_groups.empty()) fields |= P::DataField::Plddt;
    if (!config_->dssp_groups.empty()) fields |= P::DataField::Dssp;
    fields |= P::DataField::Metadata;
    return fields;
  }

private:
  std::shared_ptr<const QualityMetricsConfig> config_;
};

} // namespace

std::string normalize_group_name(std::string_view raw) {
  raw = lahuta::cli::detail::trim_view(raw);

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

void add_group_definition(std::string_view raw, bool is_plddt, std::unordered_set<std::string> &seen_names,
                          std::vector<GroupSpec> &groups) {
  raw = lahuta::cli::detail::trim_view(raw);
  if (raw.empty()) {
    throw std::runtime_error("Group definition cannot be empty.");
  }

  const auto eq_pos = raw.find('=');
  if (eq_pos == std::string_view::npos) {
    throw std::runtime_error("Group definition '" + std::string(raw) +
                             "' is missing '=' (expected name=values).");
  }

  const std::string name = normalize_group_name(raw.substr(0, eq_pos));
  if (name.empty()) {
    throw std::runtime_error("Group name '" + std::string(raw.substr(0, eq_pos)) +
                             "' is invalid (use letters, digits, '-' or '_').");
  }
  if (!seen_names.insert(name).second) {
    throw std::runtime_error("Duplicate group name '" + name + "'.");
  }

  std::vector<std::string> tokens;
  lahuta::cli::detail::split_argument_list(raw.substr(eq_pos + 1), false, tokens);
  if (tokens.empty()) {
    throw std::runtime_error("Group '" + name + "' has no values.");
  }

  std::uint32_t mask = 0;
  for (const auto &token : tokens) {
    if (is_plddt) {
      pLDDTCategory cat{};
      if (!parse_plddt_token(token, cat)) {
        throw std::runtime_error("Invalid pLDDT token '" + token + "' in group '" + name + "'.");
      }
      mask |= plddt_mask_for(cat);
    } else {
      DSSPAssignment cat{};
      if (!parse_dssp_token(token, cat)) {
        throw std::runtime_error("Invalid DSSP token '" + token + "' in group '" + name + "'.");
      }
      mask |= dssp_mask_for(cat);
    }
  }

  if (mask == 0) {
    throw std::runtime_error("Group '" + name + "' resulted in an empty mask.");
  }

  GroupSpec group;
  group.name = name;
  group.mask = mask;
  groups.push_back(std::move(group));
}

void add_default_plddt_groups(std::vector<GroupSpec> &groups) {
  GroupSpec poor;
  poor.name = "poor";
  poor.mask = plddt_mask_for(pLDDTCategory::VeryLow);
  groups.push_back(std::move(poor));

  GroupSpec low_or_poor;
  low_or_poor.name = "low_or_poor";
  low_or_poor.mask = plddt_mask_for(pLDDTCategory::Low) | plddt_mask_for(pLDDTCategory::VeryLow);
  groups.push_back(std::move(low_or_poor));
}

void add_default_dssp_groups(std::vector<GroupSpec> &groups) {
  GroupSpec coil;
  coil.name = "coil";
  coil.mask = dssp_mask_for(DSSPAssignment::Coil) | dssp_mask_for(DSSPAssignment::Turn) |
              dssp_mask_for(DSSPAssignment::Bend);
  groups.push_back(std::move(coil));

  GroupSpec helix;
  helix.name = "helix";
  helix.mask = dssp_mask_for(DSSPAssignment::AlphaHelix) | dssp_mask_for(DSSPAssignment::Helix3_10) |
               dssp_mask_for(DSSPAssignment::HelixPi) | dssp_mask_for(DSSPAssignment::PolyProlineHelix);
  groups.push_back(std::move(helix));

  GroupSpec strand;
  strand.name = "strand";
  strand.mask = dssp_mask_for(DSSPAssignment::Strand);
  groups.push_back(std::move(strand));
}

void finalize_group_keys(std::vector<GroupSpec> &plddt_groups, std::vector<GroupSpec> &dssp_groups) {
  for (auto &group : plddt_groups) {
    build_plddt_group_keys(group);
  }
  for (auto &group : dssp_groups) {
    build_dssp_group_keys(group);
  }
}

std::vector<OverlapSpec> build_overlap_specs(const std::vector<GroupSpec> &plddt_groups,
                                             const std::vector<GroupSpec> &dssp_groups) {
  std::vector<OverlapSpec> overlaps;
  overlaps.reserve(plddt_groups.size() * dssp_groups.size());
  for (std::size_t i = 0; i < plddt_groups.size(); ++i) {
    for (std::size_t j = 0; j < dssp_groups.size(); ++j) {
      OverlapSpec overlap;
      overlap.plddt_index = i;
      overlap.dssp_index  = j;
      overlap.key         = build_overlap_key(plddt_groups[i].name, dssp_groups[j].name);
      overlaps.push_back(std::move(overlap));
    }
  }
  return overlaps;
}

std::shared_ptr<P::ITask> make_quality_metrics_task(std::shared_ptr<const QualityMetricsConfig> config) {
  return std::make_shared<QualityMetricsTask>(std::move(config));
}

} // namespace lahuta::cli::quality_metrics
