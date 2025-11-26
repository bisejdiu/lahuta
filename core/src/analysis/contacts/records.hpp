#ifndef LAHUTA_ANALYSIS_CONTACTS_RECORDS_HPP
#define LAHUTA_ANALYSIS_CONTACTS_RECORDS_HPP

#include <memory>
#include <optional>
#include <string>

#include "analysis/contacts/provider.hpp"
#include "entities/contact.hpp"
#include "entities/interaction_types.hpp"
#include "topology.hpp"

namespace lahuta::analysis::contacts {

struct ContactsRecord {
  bool success;
  std::string file_path;
  std::optional<std::string> trajectory_file;
  contacts::ContactProvider provider;
  InteractionTypeSet contact_types;
  ContactSet contacts;
  std::size_t num_contacts;
  std::size_t frame_index = 0;
  std::shared_ptr<const Topology> topology;
};

} // namespace lahuta::analysis::contacts

#endif // LAHUTA_ANALYSIS_CONTACTS_RECORDS_HPP
