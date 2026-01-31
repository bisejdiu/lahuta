#include <array>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "analysis/compaction/rg_utils.hpp"
#include "analysis/compaction/shape_metrics_utils.hpp"
#include "analysis/extract/extract_tasks.hpp"
#include "logging/logging.hpp"
#include "pipeline/data/data_requirements.hpp"
#include "serialization/json.hpp"
#include "tasks/shape_metrics_task.hpp"
#include "utils/span.hpp"

namespace lahuta::cli::shape_metrics {
namespace A = lahuta::analysis;
namespace P = lahuta::pipeline;
namespace {

using lahuta::span;

class ShapeMetricsTask final : public P::ITask {
public:
  explicit ShapeMetricsTask(std::shared_ptr<const ShapeMetricsConfig> config)
      : config_(std::move(config)),
        counters_((config_ && config_->counters) ? config_->counters
                                                 : std::make_shared<ShapeMetricsCounters>()) {}

  P::TaskResult run(const std::string &item_path, P::TaskContext &ctx) override {
    const std::string *sequence              = nullptr;
    const std::vector<pLDDTCategory> *plddts = nullptr;
    const std::vector<DSSPAssignment> *dssp  = nullptr;
    const RDGeom::POINT3D_VECT *positions    = nullptr;
    std::shared_ptr<const ModelParserResult> parsed;

    auto payload = ctx.model_payload();
    if (payload) {
      if (payload->sequence && !payload->sequence->empty()) {
        sequence = payload->sequence.get();
      }
      if (payload->plddts && !payload->plddts->empty()) {
        plddts = payload->plddts.get();
      }
      if (payload->dssp && !payload->dssp->empty()) {
        dssp = payload->dssp.get();
      }
      if (payload->positions && !payload->positions->empty()) {
        positions = payload->positions.get();
      }
    } else {
      parsed = A::get_cached_model_parser_result(ctx);
      if (parsed) {
        if (!parsed->sequence.empty()) sequence = &parsed->sequence;
        if (!parsed->plddt_per_residue.empty()) plddts = &parsed->plddt_per_residue;
        if (!parsed->dssp_per_residue.empty()) dssp = &parsed->dssp_per_residue;
        if (parsed->coords_size() > 0) positions = &parsed->coords;
      }
    }

    counters_->bump_processed();

    if (!plddts) {
      Logger::get_logger()->warn("[shape-metrics] Missing pLDDT data for '{}'", item_path);
      return {};
    }

    auto &local = counters_->local();
    if (!sequence || sequence->empty()) {
      ++local.missing_sequence;
      return {};
    }
    if (!dssp || dssp->empty()) {
      ++local.missing_dssp;
      return {};
    }
    if (!positions || positions->empty()) {
      ++local.missing_positions;
      return {};
    }

    const std::size_t length = sequence->size();
    if (plddts->size() != length || dssp->size() != length) {
      ++local.length_mismatch;
      return {};
    }

    std::vector<std::size_t> atom_counts;
    try {
      atom_counts = A::atom_counts_for_sequence(*sequence);
    } catch (const std::exception &) {
      ++local.prefix_errors;
      return {};
    }

    const auto prefix = A::prefix_sums(span<const std::size_t>(atom_counts));
    if (prefix.empty()) {
      ++local.prefix_errors;
      return {};
    }
    const std::size_t raw_atom_count = prefix.back();

    if (positions->size() != raw_atom_count) {
      ++local.atom_mismatch;
      return {};
    }

    const auto coords_span = span<const RDGeom::Point3D>(*positions);
    const double raw_rg    = A::radius_of_gyration(coords_span);

    const auto trim_result = A::trim_low_confidence_tails(span<const pLDDTCategory>(*plddts),
                                                          span<const DSSPAssignment>(*dssp));

    if (trim_result.trimmed_length() == 0) {
      ++local.trimmed_empty;
      return {};
    }

    const std::size_t start_atom = prefix[trim_result.start_res];
    const std::size_t end_atom   = prefix[trim_result.end_res];
    if (end_atom <= start_atom) {
      ++local.trimmed_empty;
      return {};
    }

    const std::size_t trimmed_atom_count = end_atom - start_atom;
    const auto trimmed_coords            = span<const RDGeom::Point3D>(coords_span.data() + start_atom,
                                                            trimmed_atom_count);

    const double high_fraction_full = A::confidence_fraction(span<const pLDDTCategory>(*plddts));
    std::array<bool, A::ThresholdCount> passes{};
    for (std::size_t i = 0; i < A::ConfidenceThresholds.size(); ++i) {
      const double thr = A::ConfidenceThresholds[i];
      passes[i]        = trim_result.high_fraction >= thr;
      if (passes[i]) {
        local.add_pass(i);
      }
    }
    const double max_passed_conf = A::highest_passed_threshold(trim_result.high_fraction,
                                                               span<const double>(A::ConfidenceThresholds));

    local.total_kept_residues    += trim_result.trimmed_length();
    local.total_removed_residues += trim_result.trimmed_n + trim_result.trimmed_c;
    local.total_trimmed_n        += trim_result.trimmed_n;
    local.total_trimmed_c        += trim_result.trimmed_c;

    const double min_high_fraction = config_ ? config_->min_high_fraction : 0.80;
    if (trim_result.high_fraction < min_high_fraction) {
      counters_->bump_failed_confidence();
      return {};
    }

    double rg_trimmed      = 0.0;
    double rg_sq_trimmed   = 0.0;
    double asphericity     = 0.0;
    double acylindricity   = 0.0;
    double kappa2          = 0.0;
    double trimmed_rg_norm = 0.0;
    try {
      std::vector<double> masses(trimmed_atom_count, 1.0);
      rg_trimmed    = A::radius_of_gyration(trimmed_coords, span<const double>(masses));
      rg_sq_trimmed = A::radius_of_gyration_squared(trimmed_coords, span<const double>(masses));
      asphericity   = A::asphericity(trimmed_coords, span<const double>(masses));
      acylindricity = A::acylindricity(trimmed_coords, span<const double>(masses));
      kappa2        = A::relative_shape_anisotropy(trimmed_coords, span<const double>(masses));

      const double trimmed_rg_unweighted = A::radius_of_gyration(trimmed_coords);
      trimmed_rg_norm = A::normalized_rg(trimmed_rg_unweighted, trim_result.trimmed_length());
    } catch (const std::exception &e) {
      Logger::get_logger()->warn("[shape-metrics] Failed to compute shape metrics for '{}': {}",
                                 item_path,
                                 e.what());
      return {};
    }

    // clang-format off
    JsonBuilder json(512);
    json.key("model").value(item_path)
        .key("length").value(static_cast<std::uint64_t>(length))
        .key("trimmed_length").value(static_cast<std::uint64_t>(trim_result.trimmed_length()))
        .key("trimmed_fraction").value(length ? static_cast<double>(trim_result.trimmed_length()) / static_cast<double>(length) : 0.0)
        .key("trimmed_n").value(static_cast<std::uint64_t>(trim_result.trimmed_n))
        .key("trimmed_c").value(static_cast<std::uint64_t>(trim_result.trimmed_c))
        .key("atom_count").value(static_cast<std::uint64_t>(raw_atom_count))
        .key("trimmed_atom_count").value(static_cast<std::uint64_t>(trimmed_atom_count))
        .key("trimmed_atom_fraction")
        .value(raw_atom_count ? static_cast<double>(trimmed_atom_count) / static_cast<double>(raw_atom_count) : 0.0)
        .key("rg_raw").value(raw_rg)
        .key("rg_trimmed").value(rg_trimmed)
        .key("rg_trimmed_over_n13").value(trimmed_rg_norm)
        .key("rg_sq_trimmed").value(rg_sq_trimmed)
        .key("asphericity").value(asphericity)
        .key("acylindricity").value(acylindricity)
        .key("kappa2").value(kappa2)
        .key("plddt_high_fraction_full").value(high_fraction_full)
        .key("plddt_high_fraction_trimmed").value(trim_result.high_fraction)
        .key("max_passed_confidence").value(max_passed_conf)
        .key("passes_confidence")
        .begin_object();
    // clang-format on

    for (std::size_t i = 0; i < passes.size(); ++i) {
      json.key(A::ThresholdLabels[i]).value(passes[i]);
    }

    json.end_object();

    P::TaskResult result;
    result.ok = true;
    result.emits.push_back(P::Emission{std::string(OutputChannel), json.str()});

    counters_->bump_written();

    return result;
  }

  P::DataFieldSet data_requirements() const override {
    return P::DataFieldSet::of(
        {P::DataField::Sequence, P::DataField::Positions, P::DataField::Plddt, P::DataField::Dssp});
  }

private:
  std::shared_ptr<const ShapeMetricsConfig> config_;
  std::shared_ptr<ShapeMetricsCounters> counters_;
};

} // namespace

std::shared_ptr<P::ITask> make_shape_metrics_task(std::shared_ptr<const ShapeMetricsConfig> config) {
  return std::make_shared<ShapeMetricsTask>(std::move(config));
}

} // namespace lahuta::cli::shape_metrics
