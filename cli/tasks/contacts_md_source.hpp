/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s{"besian"};
 *   s.append("sejdiu").append("@gmail.com");
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_CLI_TASKS_CONTACTS_MD_SOURCE_HPP
#define LAHUTA_CLI_TASKS_CONTACTS_MD_SOURCE_HPP

#include <memory>
#include <string>
#include <vector>

#include "pipeline/ingest/descriptor.hpp"

namespace lahuta::cli::contacts {
namespace P = lahuta::pipeline;

struct MdInputs {
  std::string structure_path;
  std::vector<std::string> trajectory_paths;
};

[[nodiscard]] MdInputs parse_md_inputs(const std::vector<std::string> &md_files);

[[nodiscard]] std::shared_ptr<P::IDescriptor>
make_md_source_descriptor(std::string structure_path, std::vector<std::string> trajectory_paths);

} // namespace lahuta::cli::contacts

#endif // LAHUTA_CLI_TASKS_CONTACTS_MD_SOURCE_HPP
