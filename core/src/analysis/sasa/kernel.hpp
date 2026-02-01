#ifndef LAHUTA_ANALYSIS_SASA_KERNEL_HPP
#define LAHUTA_ANALYSIS_SASA_KERNEL_HPP

#include <charconv>
#include <cmath>
#include <cstdint>
#include <exception>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "analysis/compaction/rg_utils.hpp"
#include "analysis/extract/extract_tasks.hpp"
#include "analysis/sasa/protor_radii.hpp"
#include "analysis/sasa/records.hpp"
#include "analysis/sasa/sasa.hpp"
#include "compute/result.hpp"
#include "models/tables.hpp"
#include "pipeline/task/compute/context.hpp"
#include "pipeline/task/compute/parameters.hpp"
#include "pipeline/task/emission.hpp"
#include "serialization/serializer.hpp"
#include "utils/span.hpp"

namespace lahuta::analysis {
namespace P = lahuta::pipeline;
namespace C = lahuta::compute;

namespace {

using lahuta::span;

struct ThreadScratch {
  std::vector<double> radii;
  std::vector<double> per_atom;
  std::vector<std::string> labels;
};

ThreadScratch &thread_scratch() {
  static thread_local ThreadScratch scratch;
  return scratch;
}

inline void append_uint(std::string &out, std::uint64_t value) {
  char buffer[32];
  const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value);
  out.append(buffer, static_cast<std::size_t>(result.ptr - buffer));
}

inline void append_int(std::string &out, std::int64_t value) {
  char buffer[32];
  const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value);
  out.append(buffer, static_cast<std::size_t>(result.ptr - buffer));
}

[[nodiscard]] bool build_atom_radii(std::string_view sequence, std::vector<double> &out, std::string &error) {
  out.clear();
  if (sequence.empty()) {
    error = "Missing sequence data.";
    return false;
  }

  std::size_t total = 0;
  for (std::size_t idx = 0; idx < sequence.size(); ++idx) {
    const char aa = sequence[idx];
    if (!StandardAminoAcidDataTable.is_valid(aa)) {
      error = "Unsupported residue '" + std::string(1, aa) + "' at position " + std::to_string(idx);
      return false;
    }
    const auto &entry  = StandardAminoAcidDataTable[aa];
    total             += entry.size;
  }
  if (!sequence.empty()) {
    total += 1; // terminal OXT
  }
  out.reserve(total);

  for (std::size_t idx = 0; idx < sequence.size(); ++idx) {
    const char aa = sequence[idx];
    if (!StandardAminoAcidDataTable.is_valid(aa)) {
      error = "Unsupported residue '" + std::string(1, aa) + "' at position " + std::to_string(idx);
      return false;
    }
    const auto &entry = StandardAminoAcidDataTable[aa];
    if (!entry.name) {
      error = "Missing residue name for '" + std::string(1, aa) + "'";
      return false;
    }
    const std::string_view residue_name(entry.name);
    for (std::size_t atom_idx = 0; atom_idx < entry.size; ++atom_idx) {
      const char *atom_name = entry.atoms[atom_idx];
      if (!atom_name || atom_name[0] == '\0') {
        error = "Invalid atom name for residue '" + std::string(1, aa) + "'";
        return false;
      }
      const double radius = protor_radius_for_atom(residue_name, atom_name);
      if (!std::isfinite(radius) || radius <= 0.0) {
        error = "ProtOr radius missing for residue '" + std::string(residue_name) + "' atom '" + atom_name +
                "'";
        return false;
      }
      out.push_back(radius);
    }
  }

  const double oxt_radius = ProtorOxtRadius;
  if (!std::isfinite(oxt_radius) || oxt_radius <= 0.0) {
    error = "Unsupported element for OXT";
    return false;
  }
  out.push_back(oxt_radius);
  return true;
}

[[nodiscard]] bool build_atom_labels(std::string_view sequence, std::vector<std::string> &out,
                                     std::string &error) {
  out.clear();
  if (sequence.empty()) {
    error = "Missing sequence data.";
    return false;
  }

  std::size_t total = 0;
  for (std::size_t idx = 0; idx < sequence.size(); ++idx) {
    const char aa = sequence[idx];
    if (!StandardAminoAcidDataTable.is_valid(aa)) {
      error = "Unsupported residue '" + std::string(1, aa) + "' at position " + std::to_string(idx);
      return false;
    }
    const auto &entry = StandardAminoAcidDataTable[aa];
    if (!entry.name) {
      error = "Missing residue name for '" + std::string(1, aa) + "'";
      return false;
    }
    total += entry.size;
  }
  if (!sequence.empty()) {
    total += 1; // terminal OXT
  }
  out.reserve(total);

  std::size_t atom_index = 0;
  for (std::size_t idx = 0; idx < sequence.size(); ++idx) {
    const char aa        = sequence[idx];
    const auto &entry    = StandardAminoAcidDataTable[aa];
    const int residue_id = static_cast<int>(idx + 1);
    const std::string_view residue_name(entry.name);
    constexpr std::string_view chain_id = "A";

    for (std::size_t atom_idx = 0; atom_idx < entry.size; ++atom_idx) {
      const char *atom_name = entry.atoms[atom_idx];
      if (!atom_name || atom_name[0] == '\0') {
        error = "Invalid atom name for residue '" + std::string(1, aa) + "'";
        return false;
      }
      out.emplace_back();
      auto &label = out.back();
      label.reserve(32 + residue_name.size());
      append_uint(label, static_cast<std::uint64_t>(atom_index));
      label.push_back('-');
      label.append(atom_name);
      label.push_back('-');
      append_int(label, static_cast<std::int64_t>(residue_id));
      label.push_back('-');
      label.append(residue_name.data(), residue_name.size());
      label.push_back('-');
      label.append(chain_id.data(), chain_id.size());
      ++atom_index;
    }
  }

  const std::size_t last_idx = sequence.size() - 1;
  const auto &last_entry     = StandardAminoAcidDataTable[sequence[last_idx]];
  const std::string_view last_residue_name(last_entry.name);
  const int residue_id                = static_cast<int>(sequence.size());
  constexpr std::string_view chain_id = "A";
  out.emplace_back();
  auto &label = out.back();
  label.reserve(32 + last_residue_name.size());
  append_uint(label, static_cast<std::uint64_t>(atom_index));
  label.push_back('-');
  label.append("OXT");
  label.push_back('-');
  append_int(label, static_cast<std::int64_t>(residue_id));
  label.push_back('-');
  label.append(last_residue_name.data(), last_residue_name.size());
  label.push_back('-');
  label.append(chain_id.data(), chain_id.size());
  return true;
}

} // namespace

struct SasaSrKernel {
  static C::ComputationResult execute(C::DataContext<P::PipelineContext, C::Mut::ReadWrite> &context,
                                      const P::SasaSrParams &p) {
    auto &data = context.data();
    std::string_view sequence;
    const RDGeom::POINT3D_VECT *positions = nullptr;
    std::shared_ptr<const ModelParserResult> parsed;

    auto model_payload = data.ctx ? data.ctx->model_payload() : nullptr;
    if (model_payload) {
      if (model_payload->sequence && !model_payload->sequence->empty()) {
        sequence = *model_payload->sequence;
      } else if (model_payload->sequence_view && !model_payload->sequence_view->data.empty()) {
        sequence = model_payload->sequence_view->data;
      }
      if (model_payload->positions && !model_payload->positions->empty()) {
        positions = model_payload->positions.get();
      }
    } else if (data.ctx) {
      parsed = get_cached_model_parser_result(*data.ctx);
      if (parsed) {
        if (!parsed->sequence.empty()) sequence = parsed->sequence;
        if (parsed->coords_size() > 0) positions = &parsed->coords;
      }
    }

    auto counters = p.counters;
    if (!counters) {
      static auto fallback = std::make_shared<SasaSrCounters>();
      counters             = fallback;
    }

    const auto processed = counters->processed.fetch_add(1, std::memory_order_relaxed) + 1;

    if (sequence.empty()) {
      counters->missing_sequence.fetch_add(1, std::memory_order_relaxed);
      return C::ComputationResult(
          C::ComputationError("SASA-SR missing sequence for '" + data.item_path + "'"));
    }
    if (!positions || positions->empty()) {
      counters->missing_positions.fetch_add(1, std::memory_order_relaxed);
      return C::ComputationResult(
          C::ComputationError("SASA-SR missing positions for '" + data.item_path + "'"));
    }

    std::size_t expected_atoms = 0;
    try {
      const auto counts = atom_counts_for_sequence(sequence);
      for (const auto count : counts) {
        expected_atoms += count;
      }
    } catch (const std::exception &e) {
      counters->invalid_residue.fetch_add(1, std::memory_order_relaxed);
      return C::ComputationResult(C::ComputationError("SASA-SR invalid sequence for '" + data.item_path +
                                                      "': " + std::string(e.what())));
    }

    if (positions->size() != expected_atoms) {
      counters->atom_mismatch.fetch_add(1, std::memory_order_relaxed);
      return C::ComputationResult(C::ComputationError("SASA-SR atom count mismatch for '" + data.item_path +
                                                      "': positions=" + std::to_string(positions->size()) +
                                                      " expected=" + std::to_string(expected_atoms)));
    }

    auto &scratch = thread_scratch();
    std::string error;
    if (!build_atom_radii(sequence, scratch.radii, error)) {
      counters->invalid_residue.fetch_add(1, std::memory_order_relaxed);
      return C::ComputationResult(
          C::ComputationError("SASA-SR radii mapping failed for '" + data.item_path + "': " + error));
    }
    if (!build_atom_labels(sequence, scratch.labels, error)) {
      counters->invalid_residue.fetch_add(1, std::memory_order_relaxed);
      return C::ComputationResult(
          C::ComputationError("SASA-SR label mapping failed for '" + data.item_path + "': " + error));
    }
    if (scratch.radii.size() != positions->size()) {
      counters->atom_mismatch.fetch_add(1, std::memory_order_relaxed);
      return C::ComputationResult(C::ComputationError("SASA-SR radii mismatch for '" + data.item_path +
                                                      "': radii=" + std::to_string(scratch.radii.size()) +
                                                      " positions=" + std::to_string(positions->size())));
    }
    if (scratch.labels.size() != positions->size()) {
      counters->atom_mismatch.fetch_add(1, std::memory_order_relaxed);
      return C::ComputationResult(C::ComputationError("SASA-SR label mismatch for '" + data.item_path +
                                                      "': labels=" + std::to_string(scratch.labels.size()) +
                                                      " positions=" + std::to_string(positions->size())));
    }

    SasaParams local_params = p.params;

    const auto coords_span = span<const RDGeom::Point3D>(*positions);
    AtomView atoms;
    atoms.coords = coords_span;
    atoms.radii  = span<const double>(scratch.radii);

    auto result = compute_sasa(atoms, local_params);

    SasaSrRecord record;
    record.model_path    = data.item_path;
    record.labels        = std::move(scratch.labels);
    record.per_atom      = std::move(result.per_atom);
    record.total         = result.total;
    record.include_total = p.include_total;

    auto payload = serialization::Serializer<fmt::json, SasaSrRecord>::serialize(record);

    scratch.labels   = std::move(record.labels);
    scratch.per_atom = std::move(record.per_atom);

    P::EmissionList out;
    out.push_back(P::Emission{p.channel, std::move(payload)});

    counters->written.fetch_add(1, std::memory_order_relaxed);

    return C::ComputationResult(std::move(out));
  }
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SASA_KERNEL_HPP
