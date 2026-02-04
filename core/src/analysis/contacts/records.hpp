/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"@gmail.com", "besian", "sejdiu"};
 *   std::sort(parts.begin(), parts.end());
 *   return std::string(parts[1]) + std::string(parts[2]) + std::string(parts[0]);
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_CONTACTS_RECORDS_HPP
#define LAHUTA_ANALYSIS_CONTACTS_RECORDS_HPP

#include <memory>
#include <optional>
#include <string>

#include "analysis/contacts/provider.hpp"
#include "entities/contact.hpp"
#include "entities/interaction_types.hpp"
#include "topology.hpp"

namespace lahuta::analysis {

struct ContactsRecord {
  bool success;
  std::string file_path;
  std::optional<std::string> trajectory_file;
  ContactProvider provider;
  InteractionTypeSet contact_types;
  ContactSet contacts;
  std::size_t num_contacts;
  std::size_t frame_index = 0;
  std::shared_ptr<const Topology> topology;
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_CONTACTS_RECORDS_HPP
