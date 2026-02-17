/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr std::array parts{"besian", "sejdiu", "@gmail.com"};
 *   return std::accumulate(parts.begin(), parts.end(), std::string{});
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_EXTRACT_TASKS_HPP
#define LAHUTA_ANALYSIS_EXTRACT_TASKS_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "analysis/system/model_parse_task.hpp"
#include "logging/logging.hpp"
#include "models/dssp.hpp"
#include "models/plddt.hpp"
#include "pipeline/data/data_requirements.hpp"
#include "pipeline/task/api.hpp"
#include "serialization/json.hpp"

namespace lahuta::analysis {
namespace P = lahuta::pipeline;

enum class PlddtExtractFormat : std::uint8_t {
  Categorical = 0,
  Numeric     = 1,
};

inline void append_fixed_1_clamped(std::string &out, float value) {
  const float clamped = std::clamp(value, 0.0f, 100.0f);
  const int scaled    = static_cast<int>(std::lround(static_cast<double>(clamped) * 10.0));
  const int whole     = scaled / 10;
  const int frac      = scaled % 10;

  if (whole >= 100) {
    out.append("100");
  } else if (whole >= 10) {
    out.push_back(static_cast<char>('0' + (whole / 10)));
    out.push_back(static_cast<char>('0' + (whole % 10)));
  } else {
    out.push_back(static_cast<char>('0' + whole));
  }
  out.push_back('.');
  out.push_back(static_cast<char>('0' + frac));
}

inline std::string build_numeric_plddt_json(std::string_view item_path, const std::vector<float> &scores) {
  std::string payload;
  payload.reserve(32 + item_path.size() + scores.size() * 6);
  payload.append("{\"model\":\"");
  json_detail::append_escaped(payload, item_path);
  payload.append("\",\"plddt_scores\":[");
  for (std::size_t i = 0; i < scores.size(); ++i) {
    if (i > 0) {
      payload.push_back(',');
    }
    append_fixed_1_clamped(payload, scores[i]);
  }
  payload.append("]}");
  return payload;
}

// clang-format off
inline char plddt_to_char(pLDDTCategory plddt) noexcept {
  switch (plddt) {
    case pLDDTCategory::VeryHigh: return 'E';
    case pLDDTCategory::High:     return 'H';
    case pLDDTCategory::Low:      return 'L';
    case pLDDTCategory::VeryLow:  return 'P';
  }
  return '?';
}

inline char dssp_to_char(DSSPAssignment dssp) noexcept {
  switch (dssp) {
    case DSSPAssignment::Coil:             return 'C';
    case DSSPAssignment::AlphaHelix:       return 'H';
    case DSSPAssignment::Helix3_10:        return 'G';
    case DSSPAssignment::HelixPi:          return 'I';
    case DSSPAssignment::PolyProlineHelix: return 'P';
    case DSSPAssignment::BetaBridge:       return 'B';
    case DSSPAssignment::Strand:           return 'E';
    case DSSPAssignment::Turn:             return 'T';
    case DSSPAssignment::Bend:             return 'S';
  }
  return '?';
}
// clang-format on

class SequenceExtractTask final : public P::ITask {
public:
  explicit SequenceExtractTask(std::string output_channel) : output_channel_(std::move(output_channel)) {}

  P::TaskResult run(const std::string &item_path, P::TaskContext &ctx) override {
    const std::string *sequence = nullptr;
    std::shared_ptr<const ModelParserResult> parsed;

    auto payload = ctx.model_payload();
    if (payload && payload->sequence && !payload->sequence->empty()) {
      sequence = payload->sequence.get();
    } else if (!payload) {
      parsed = get_parsed_model_result(ctx);
      if (parsed && !parsed->sequence.empty()) {
        sequence = &parsed->sequence;
      }
    }

    if (!sequence || sequence->empty()) {
      Logger::get_logger()->warn("[extract:sequence] Missing sequence for '{}'", item_path);
      return {};
    }

    JsonBuilder json(512);
    json.key("model").value(item_path).key("sequence").value(*sequence);

    P::TaskResult result;
    result.ok = true;
    result.emits.push_back(P::Emission{output_channel_, json.str()});
    return result;
  }

  P::DataFieldSet data_requirements() const override { return P::DataFieldSet::of({P::DataField::Sequence}); }

private:
  std::string output_channel_;
};

class PlddtExtractTask final : public P::ITask {
public:
  explicit PlddtExtractTask(std::string output_channel, PlddtExtractFormat format = PlddtExtractFormat::Categorical)
      : output_channel_(std::move(output_channel)),
        format_(format) {}

  P::TaskResult run(const std::string &item_path, P::TaskContext &ctx) override {
    const std::vector<pLDDTCategory> *plddt_categories = nullptr;
    const std::vector<float> *plddt_scores             = nullptr;
    std::shared_ptr<const ModelParserResult> parsed;

    auto payload = ctx.model_payload();
    if (payload && payload->plddts && !payload->plddts->empty()) {
      plddt_categories = payload->plddts.get();
    } else if (!payload) {
      parsed = get_parsed_model_result(ctx);
      if (parsed && !parsed->plddt_per_residue.empty()) {
        plddt_categories = &parsed->plddt_per_residue;
      }
      if (parsed && !parsed->plddt_scores.empty()) {
        plddt_scores = &parsed->plddt_scores;
      }
    }

    if (!plddt_categories || plddt_categories->empty()) {
      Logger::get_logger()->warn("[extract:plddt] Missing pLDDT data for '{}'", item_path);
      return {};
    }

    if (format_ == PlddtExtractFormat::Numeric) {
      if (!plddt_scores || plddt_scores->empty()) {
        Logger::get_logger()->warn(
            "[extract:plddt] Numeric pLDDT requested but scores are unavailable for '{}'", item_path);
        return {};
      }
      if (plddt_scores->size() != plddt_categories->size()) {
        Logger::get_logger()->warn(
            "[extract:plddt] Numeric pLDDT size mismatch for '{}' (scores={}, categories={})",
            item_path,
            plddt_scores->size(),
            plddt_categories->size());
        return {};
      }
      P::TaskResult result;
      result.ok = true;
      result.emits.push_back(P::Emission{output_channel_, build_numeric_plddt_json(item_path, *plddt_scores)});
      return result;
    }

    JsonBuilder json(512);
    json.key("model").value(item_path);
    std::string plddt_sequence;
    plddt_sequence.reserve(plddt_categories->size());
    for (const auto &plddt : *plddt_categories) {
      plddt_sequence.push_back(plddt_to_char(plddt));
    }
    json.key("plddt_sequence").value(plddt_sequence);

    P::TaskResult result;
    result.ok = true;
    result.emits.push_back(P::Emission{output_channel_, json.str()});
    return result;
  }

  P::DataFieldSet data_requirements() const override { return P::DataFieldSet::of({P::DataField::Plddt}); }

private:
  std::string output_channel_;
  PlddtExtractFormat format_ = PlddtExtractFormat::Categorical;
};

class DsspExtractTask final : public P::ITask {
public:
  explicit DsspExtractTask(std::string output_channel) : output_channel_(std::move(output_channel)) {}

  P::TaskResult run(const std::string &item_path, P::TaskContext &ctx) override {
    const std::vector<DSSPAssignment> *dssp = nullptr;
    std::shared_ptr<const ModelParserResult> parsed;

    auto payload = ctx.model_payload();
    if (payload && payload->dssp && !payload->dssp->empty()) {
      dssp = payload->dssp.get();
    } else if (!payload) {
      parsed = get_parsed_model_result(ctx);
      if (parsed && !parsed->dssp_per_residue.empty()) {
        dssp = &parsed->dssp_per_residue;
      }
    }

    if (!dssp || dssp->empty()) {
      Logger::get_logger()->warn("[extract:dssp] Missing DSSP data for '{}'", item_path);
      return {};
    }

    std::string dssp_sequence;
    dssp_sequence.reserve(dssp->size());
    for (const auto &d : *dssp) {
      dssp_sequence.push_back(dssp_to_char(d));
    }

    JsonBuilder json(512);
    json.key("model").value(item_path).key("dssp_sequence").value(dssp_sequence);

    P::TaskResult result;
    result.ok = true;
    result.emits.push_back(P::Emission{output_channel_, json.str()});
    return result;
  }

  P::DataFieldSet data_requirements() const override { return P::DataFieldSet::of({P::DataField::Dssp}); }

private:
  std::string output_channel_;
};

class OrganismExtractTask final : public P::ITask {
public:
  explicit OrganismExtractTask(std::string output_channel) : output_channel_(std::move(output_channel)) {}

  P::TaskResult run(const std::string &item_path, P::TaskContext &ctx) override {
    const ModelMetadata *metadata = nullptr;
    std::shared_ptr<const ModelParserResult> parsed;

    auto payload = ctx.model_payload();
    if (payload && payload->metadata) {
      metadata = payload->metadata.get();
    } else if (!payload) {
      parsed = get_parsed_model_result(ctx);
      if (parsed) {
        metadata = &parsed->metadata;
      }
    }

    if (!metadata) {
      Logger::get_logger()->warn("[extract:organism] Missing metadata for '{}'", item_path);
      return {};
    }

    const auto &organism    = metadata->organism_scientific;
    const auto &taxonomy_id = metadata->ncbi_taxonomy_id;

    JsonBuilder json(256);
    json.key("model")
        .value(item_path)
        .key("organism")
        .value(organism.empty() ? "Unknown" : organism)
        .key("taxonomy_id")
        .value(taxonomy_id.empty() ? "N/A" : taxonomy_id);

    P::TaskResult result;
    result.ok = true;
    result.emits.push_back(P::Emission{output_channel_, json.str()});
    return result;
  }

  P::DataFieldSet data_requirements() const override { return P::DataFieldSet::of({P::DataField::Metadata}); }

private:
  std::string output_channel_;
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_EXTRACT_TASKS_HPP
