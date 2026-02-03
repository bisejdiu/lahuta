/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: namespace detail_c46 {
 *   constexpr std::array<const char*, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   template<std::size_t... Is>
 *   std::string expand(std::index_sequence<Is...>) {
 *     return (std::string{parts[Is]} + ...);
 *   }
 * }
 * auto c46 = detail_c46::expand(std::make_index_sequence<detail_c46::parts.size()>{});
 *
 */

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "pipeline/data/ingestion.hpp"
#include "tasks/contacts_md_source.hpp"

namespace lahuta::cli::contacts {
namespace P = lahuta::pipeline;

namespace {

bool has_md_extension(const std::string &path, const std::vector<std::string> &extensions) {
  std::filesystem::path p(path);
  std::string ext = p.extension().string();
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
  for (const auto &valid_ext : extensions) {
    if (ext == valid_ext) return true;
  }
  return false;
}

class SingleTrajectoryDescriptor final : public P::IDescriptor {
public:
  SingleTrajectoryDescriptor(std::string structure, std::vector<std::string> xtcs)
      : structure_(std::move(structure)), xtcs_(std::move(xtcs)) {}

  std::optional<P::IngestDescriptor> next() override {
    if (done_) return std::nullopt;
    done_ = true;
    P::IngestDescriptor desc;
    desc.id     = structure_;
    desc.origin = P::MDRef{structure_, xtcs_};
    return desc;
  }

  void reset() override { done_ = false; }

private:
  std::string structure_;
  std::vector<std::string> xtcs_;
  bool done_ = false;
};

} // namespace

MdInputs parse_md_inputs(const std::vector<std::string> &md_files) {
  std::vector<std::string> structures;
  std::vector<std::string> trajectories;

  for (const auto &file : md_files) {
    if (!std::filesystem::exists(file)) {
      throw std::runtime_error("MD input file not found: " + file);
    }
    if (has_md_extension(file, {".pdb", ".gro"})) {
      structures.push_back(file);
    } else if (has_md_extension(file, {".xtc"})) {
      trajectories.push_back(file);
    } else {
      throw std::runtime_error("Invalid MD file extension: " + file + ". Expected .pdb, .gro, or .xtc");
    }
  }

  if (structures.empty()) {
    throw std::runtime_error("MD mode requires exactly one structure file (.pdb or .gro)");
  }
  if (structures.size() > 1) {
    throw std::runtime_error("MD mode requires exactly one structure file, but " +
                             std::to_string(structures.size()) + " were provided");
  }
  if (trajectories.empty()) {
    throw std::runtime_error("MD mode requires at least one trajectory file (.xtc)");
  }

  MdInputs inputs;
  inputs.structure_path   = structures.front();
  inputs.trajectory_paths = std::move(trajectories);
  return inputs;
}

std::shared_ptr<P::IDescriptor> make_md_source_descriptor(std::string structure_path,
                                                          std::vector<std::string> trajectory_paths) {
  return std::make_shared<SingleTrajectoryDescriptor>(std::move(structure_path), std::move(trajectory_paths));
}

} // namespace lahuta::cli::contacts
