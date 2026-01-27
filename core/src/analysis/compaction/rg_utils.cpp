#include <cmath>
#include <stdexcept>
#include <string>

#include "analysis/compaction/rg_utils.hpp"

namespace lahuta::analysis {
namespace {

constexpr bool is_plddt_high(pLDDTCategory cat) {
  return cat == pLDDTCategory::VeryHigh || cat == pLDDTCategory::High;
}

constexpr bool is_plddt_low(pLDDTCategory cat) {
  return cat == pLDDTCategory::Low || cat == pLDDTCategory::VeryLow;
}

constexpr bool is_dssp_coil(DSSPAssignment dssp) {
  return dssp == DSSPAssignment::Coil || dssp == DSSPAssignment::Turn || dssp == DSSPAssignment::Bend;
}

// clang-format off
std::size_t atom_count_for_residue(char aa) {
  switch (aa) {
    case 'A': return 5;
    case 'R': return 11;
    case 'N': return 8;
    case 'D': return 8;
    case 'C': return 6;
    case 'Q': return 9;
    case 'E': return 9;
    case 'G': return 4;
    case 'H': return 10;
    case 'I': return 8;
    case 'L': return 8;
    case 'K': return 9;
    case 'M': return 8;
    case 'F': return 11;
    case 'P': return 7;
    case 'S': return 6;
    case 'T': return 7;
    case 'W': return 14;
    case 'Y': return 12;
    case 'V': return 7;
    default: break;
  }
  return 0;
}
// clang-format on

} // namespace

double confidence_fraction(span<const pLDDTCategory> plddt) {
  if (plddt.empty()) return 0.0;
  std::size_t high = 0;
  for (const auto cat : plddt) {
    if (is_plddt_high(cat)) ++high;
  }
  return static_cast<double>(high) / static_cast<double>(plddt.size());
}

TrimResult trim_low_confidence_tails(span<const pLDDTCategory> plddt, span<const DSSPAssignment> dssp) {
  if (plddt.size() != dssp.size()) {
    throw std::runtime_error("pLDDT and DSSP length mismatch during trimming");
  }

  std::size_t start = 0;
  std::size_t end   = plddt.size();

  while (start < end && is_dssp_coil(dssp[start]) && is_plddt_low(plddt[start])) {
    ++start;
  }

  while (end > start && is_dssp_coil(dssp[end - 1]) && is_plddt_low(plddt[end - 1])) {
    --end;
  }

  const auto trimmed = span<const pLDDTCategory>(plddt.data() + start, end - start);
  return TrimResult{
      start,
      end,
      start,
      plddt.size() - end,
      confidence_fraction(trimmed),
  };
}

std::vector<std::size_t> atom_counts_for_sequence(std::string_view sequence) {
  std::vector<std::size_t> counts;
  counts.reserve(sequence.size());
  for (std::size_t idx = 0; idx < sequence.size(); ++idx) {
    const char aa    = sequence[idx];
    const auto count = atom_count_for_residue(aa);
    if (count == 0) {
      std::string message = "Unsupported residue '";
      message.push_back(aa);
      message += "' at position ";
      message += std::to_string(idx);
      throw std::runtime_error(message);
    }
    counts.push_back(count);
  }
  if (!counts.empty()) {
    counts.back() += 1; // terminal OXT
  }
  return counts;
}

std::vector<std::size_t> prefix_sums(span<const std::size_t> values) {
  std::vector<std::size_t> prefix;
  prefix.reserve(values.size() + 1);
  prefix.push_back(0);
  std::size_t total = 0;
  for (const auto value : values) {
    total += value;
    prefix.push_back(total);
  }
  return prefix;
}

double radius_of_gyration(span<const RDGeom::Point3D> coords) {
  if (coords.empty()) return 0.0;

  double cx = 0.0;
  double cy = 0.0;
  double cz = 0.0;
  for (const auto &pt : coords) {
    cx += pt.x;
    cy += pt.y;
    cz += pt.z;
  }
  const double inv_n = 1.0 / static_cast<double>(coords.size());

  cx *= inv_n;
  cy *= inv_n;
  cz *= inv_n;

  double sum_sq = 0.0;
  for (const auto &pt : coords) {
    const double dx  = pt.x - cx;
    const double dy  = pt.y - cy;
    const double dz  = pt.z - cz;
    sum_sq          += dx * dx + dy * dy + dz * dz;
  }
  const double rg_sq = sum_sq * inv_n;
  return std::sqrt(rg_sq);
}

double radius_of_gyration_weighted(span<const RDGeom::Point3D> coords, span<const double> weights) {
  if (coords.empty()) return 0.0;
  if (weights.size() != coords.size()) {
    throw std::runtime_error("Expected weights to match coordinate count in radius_of_gyration_weighted");
  }

  double w_sum = 0.0;
  double cx    = 0.0;
  double cy    = 0.0;
  double cz    = 0.0;
  for (std::size_t i = 0; i < coords.size(); ++i) {
    const double w = weights[i];

    w_sum += w;
    cx    += w * coords[i].x;
    cy    += w * coords[i].y;
    cz    += w * coords[i].z;
  }
  if (w_sum <= 0.0) return 0.0;

  cx /= w_sum;
  cy /= w_sum;
  cz /= w_sum;

  double sum_sq = 0.0;
  for (std::size_t i = 0; i < coords.size(); ++i) {
    const double dx  = coords[i].x - cx;
    const double dy  = coords[i].y - cy;
    const double dz  = coords[i].z - cz;
    sum_sq          += weights[i] * (dx * dx + dy * dy + dz * dz);
  }
  return std::sqrt(sum_sq / w_sum);
}

double normalized_rg(double rg, std::size_t residue_count) {
  if (residue_count == 0) return 0.0;
  return rg / std::cbrt(static_cast<double>(residue_count));
}

double highest_passed_threshold(double high_fraction, span<const double> thresholds) {
  double best = 0.0;
  for (const auto thr : thresholds) {
    if (high_fraction >= thr && thr > best) {
      best = thr;
    }
  }
  return best;
}

} // namespace lahuta::analysis
