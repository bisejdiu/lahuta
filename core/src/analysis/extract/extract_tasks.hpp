#ifndef LAHUTA_ANALYSIS_EXTRACT_TASKS_HPP
#define LAHUTA_ANALYSIS_EXTRACT_TASKS_HPP

#include <memory>
#include <string>
#include <vector>

#include "analysis/system/model_loader.hpp"
#include "logging/logging.hpp"
#include "models/dssp.hpp"
#include "models/plddt.hpp"
#include "pipeline/data/data_requirements.hpp"
#include "pipeline/task/api.hpp"
#include "serialization/json.hpp"

namespace lahuta::analysis {
namespace P = lahuta::pipeline;

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
    case DSSPAssignment::Strand:           return 'E';
    case DSSPAssignment::Turn:             return 'T';
    case DSSPAssignment::Bend:             return 'S';
  }
  return '?';
}

// clang-format on
inline constexpr const char *CTX_PARSED_MODEL_KEY = "lahuta.parsed_model";

inline std::shared_ptr<const ModelParserResult> get_cached_model_parser_result(const P::TaskContext &ctx) {
  return ctx.get_object<ModelParserResult>(CTX_PARSED_MODEL_KEY);
}

class ModelParseTask final : public P::ITask {
public:
  P::TaskResult run(const std::string &item_path, P::TaskContext &ctx) override {
    if (ctx.model_payload()) return {};
    if (get_cached_model_parser_result(ctx)) return {};

    try {
      auto parsed = load_model_parser_result(item_path);
      ctx.set_object<ModelParserResult>(CTX_PARSED_MODEL_KEY,
                                        std::make_shared<ModelParserResult>(std::move(parsed)));
    } catch (const std::exception &e) {
      Logger::get_logger()->warn("[extract:parse] Failed to parse '{}': {}", item_path, e.what());
    }
    return {};
  }
};

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
      parsed = get_cached_model_parser_result(ctx);
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
  explicit PlddtExtractTask(std::string output_channel) : output_channel_(std::move(output_channel)) {}

  P::TaskResult run(const std::string &item_path, P::TaskContext &ctx) override {
    const std::vector<pLDDTCategory> *plddts = nullptr;
    std::shared_ptr<const ModelParserResult> parsed;

    auto payload = ctx.model_payload();
    if (payload && payload->plddts && !payload->plddts->empty()) {
      plddts = payload->plddts.get();
    } else if (!payload) {
      parsed = get_cached_model_parser_result(ctx);
      if (parsed && !parsed->plddt_per_residue.empty()) {
        plddts = &parsed->plddt_per_residue;
      }
    }

    if (!plddts || plddts->empty()) {
      Logger::get_logger()->warn("[extract:plddt] Missing pLDDT data for '{}'", item_path);
      return {};
    }

    std::string plddt_sequence;
    plddt_sequence.reserve(plddts->size());
    for (const auto &plddt : *plddts) {
      plddt_sequence.push_back(plddt_to_char(plddt));
    }

    JsonBuilder json(512);
    json.key("model").value(item_path).key("plddt_sequence").value(plddt_sequence);

    P::TaskResult result;
    result.ok = true;
    result.emits.push_back(P::Emission{output_channel_, json.str()});
    return result;
  }

  P::DataFieldSet data_requirements() const override { return P::DataFieldSet::of({P::DataField::Plddt}); }

private:
  std::string output_channel_;
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
      parsed = get_cached_model_parser_result(ctx);
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
      parsed = get_cached_model_parser_result(ctx);
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
