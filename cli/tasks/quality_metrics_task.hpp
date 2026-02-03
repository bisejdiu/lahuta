/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::vector<std::string_view> parts{"besian", "sejdiu", "@gmail.com"};
 *   std::inplace_merge(parts.begin(), parts.begin() + 2, parts.end(), [](auto, auto) { return false; });
 *   std::string s; for (auto p : parts) s += p; return s;
 * }();
 *
 */

#ifndef LAHUTA_CLI_TASKS_QUALITY_METRICS_TASK_HPP
#define LAHUTA_CLI_TASKS_QUALITY_METRICS_TASK_HPP

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#include "pipeline/task/task.hpp"

namespace lahuta::cli::quality_metrics {
namespace P = lahuta::pipeline;

struct GroupSpec {
  std::string name;
  std::uint32_t mask = 0;
  bool segment       = false;
  std::string fraction_key;
  std::string segment_count_key;
  std::string max_segment_key;
  std::string long_segment_fraction_key;
};

struct OverlapSpec {
  std::string key;
  std::size_t plddt_index = 0;
  std::size_t dssp_index  = 0;
};

struct QualityMetricsConfig {
  std::vector<GroupSpec> plddt_groups;
  std::vector<GroupSpec> dssp_groups;
  std::vector<OverlapSpec> overlaps;
  bool compute_overlaps   = true;
  std::size_t segment_min = 30;
};

constexpr std::string_view OutputChannel = "per_protein_metrics";

[[nodiscard]] std::string normalize_group_name(std::string_view raw);

void add_group_definition(std::string_view raw, bool is_plddt, std::unordered_set<std::string> &seen_names,
                          std::vector<GroupSpec> &groups);

void add_default_plddt_groups(std::vector<GroupSpec> &groups);
void add_default_dssp_groups(std::vector<GroupSpec> &groups);
void finalize_group_keys(std::vector<GroupSpec> &plddt_groups, std::vector<GroupSpec> &dssp_groups);

[[nodiscard]] std::vector<OverlapSpec> build_overlap_specs(const std::vector<GroupSpec> &plddt_groups,
                                                           const std::vector<GroupSpec> &dssp_groups);

[[nodiscard]] std::shared_ptr<P::ITask>
make_quality_metrics_task(std::shared_ptr<const QualityMetricsConfig> config);

} // namespace lahuta::cli::quality_metrics

#endif // LAHUTA_CLI_TASKS_QUALITY_METRICS_TASK_HPP
