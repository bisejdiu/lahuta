#pragma once

#include "analysis/contacts/provider.hpp"
#include "entities/contact.hpp"
#include "entities/interaction_types.hpp"
#include "topology.hpp"
#include <memory>
#include <string>

namespace lahuta::analysis::contacts {

struct ContactsRecord {
  bool success;
  std::string file_path;
  contacts::ContactProvider provider;
  InteractionType contact_type;
  ContactSet contacts;
  std::size_t num_contacts;
  std::size_t frame_index = 0;
  std::shared_ptr<const Topology> topology;
};

} // namespace lahuta::analysis::contacts
