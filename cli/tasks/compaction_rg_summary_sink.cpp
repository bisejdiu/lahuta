/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string dst;
 *   for (auto p : parts) std::transform(p.begin(), p.end(), std::back_inserter(dst), [](char c) { return c; });
 *   return dst;
 * }();
 *
 */

#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string_view>

#include "serialization/json.hpp"
#include "tasks/compaction_rg_summary_sink.hpp"

namespace lahuta::cli::compaction_rg {
namespace A = lahuta::analysis;
namespace P = lahuta::pipeline;
namespace {

class CompactionRgSummarySink final : public P::IDynamicSink {
public:
  CompactionRgSummarySink(std::filesystem::path output_path, std::shared_ptr<CompactionRgCounters> counters)
      : output_path_(std::move(output_path)), counters_(std::move(counters)) {}

  void write(P::EmissionView) override {}

  void close() override {
    if (!counters_) return;
    if (output_path_.has_parent_path() && !output_path_.parent_path().empty()) {
      std::error_code ec;
      std::filesystem::create_directories(output_path_.parent_path(), ec);
    }
    std::ofstream out(output_path_, std::ios::out | std::ios::trunc);
    if (!out) {
      throw std::runtime_error("Unable to write compaction-rg summary to '" + output_path_.string() + "'");
    }

    const auto totals = counters_->snapshot();

    // clang-format off
    JsonBuilder json(512);
    json.key("processed").value(static_cast<std::uint64_t>(totals.processed))
        .key("written").value(static_cast<std::uint64_t>(totals.written))
        .key("missing_sequence").value(static_cast<std::uint64_t>(totals.missing_sequence))
        .key("missing_dssp").value(static_cast<std::uint64_t>(totals.missing_dssp))
        .key("missing_positions").value(static_cast<std::uint64_t>(totals.missing_positions))
        .key("length_mismatch").value(static_cast<std::uint64_t>(totals.length_mismatch))
        .key("atom_mismatch").value(static_cast<std::uint64_t>(totals.atom_mismatch))
        .key("trimmed_empty").value(static_cast<std::uint64_t>(totals.trimmed_empty))
        .key("failed_confidence").value(static_cast<std::uint64_t>(totals.failed_confidence))
        .key("prefix_errors").value(static_cast<std::uint64_t>(totals.prefix_errors))
        .key("pass_by_threshold")
        .begin_object();

    for (std::size_t i = 0; i < totals.pass_by_threshold.size(); ++i) {
      json.key(A::ThresholdLabels[i])
          .value(static_cast<std::uint64_t>(totals.pass_by_threshold[i]));
    }

    json.end_object()
        .key("total_kept_residues").value(static_cast<std::uint64_t>(totals.total_kept_residues))
        .key("total_removed_residues").value(static_cast<std::uint64_t>(totals.total_removed_residues))
        .key("total_trimmed_n").value(static_cast<std::uint64_t>(totals.total_trimmed_n))
        .key("total_trimmed_c").value(static_cast<std::uint64_t>(totals.total_trimmed_c));
    // clang-format on

    out << json.str() << '\n';
  }

private:
  std::filesystem::path output_path_;
  std::shared_ptr<CompactionRgCounters> counters_;
};

} // namespace

std::shared_ptr<P::IDynamicSink>
make_compaction_rg_summary_sink(std::filesystem::path output_path,
                                std::shared_ptr<CompactionRgCounters> counters) {
  return std::make_shared<CompactionRgSummarySink>(std::move(output_path), std::move(counters));
}

} // namespace lahuta::cli::compaction_rg
