#ifndef LAHUTA_PIPELINE_TASKS_CONTACTS_TASK_HPP
#define LAHUTA_PIPELINE_TASKS_CONTACTS_TASK_HPP

#include "contacts/arpeggio/provider.hpp"
#include "contacts/engine.hpp"
#include "contacts/molstar/provider.hpp"
#include "entities/contact.hpp"
#include "entities/interaction_types.hpp"
#include "lahuta.hpp"
#include "logging.hpp"
#include "topology.hpp"
#include <memory>
#include <string>
#include <string_view>

// clang-format off
namespace lahuta::tasks {

enum class ContactProvider { MolStar, Arpeggio };

class ContactsTask {
public:
  using input_type = std::string_view;

  struct result_type {
    bool success;
    std::string file_path;
    ContactProvider provider;
    InteractionType contact_type;  // we'll use none for all contacts
    ContactSet contacts;
    size_t num_contacts;
    std::unique_ptr<const Topology> topology;
  };

  explicit ContactsTask(ContactProvider provider = ContactProvider::MolStar, InteractionType contact_type = InteractionType::None)
    : provider_(provider), contact_type_(contact_type) {}

  result_type operator()(std::string_view file_path) const {
    result_type result;
    result.file_path    = std::string(file_path);
    result.provider     = provider_;
    result.contact_type = contact_type_;
    result.success      = false;
    result.num_contacts = 0;
    result.topology     = nullptr;

    try {
      Luni luni(result.file_path);

      const auto computer_type = (provider_ == ContactProvider::Arpeggio)
        ? ContactComputerType::Arpeggio
        : ContactComputerType::Molstar;

      luni.set_atom_typing_method(computer_type);

      if (!luni.build_topology()) {
        Logger::get_logger()->error("Failed to build topology for file: {}", file_path);
        return result;
      }

      // Luni will become invalid when we return so we take ownership of the topology
      result.topology = luni.release_topology();

      if (provider_ == ContactProvider::Arpeggio) {
        InteractionEngine<ArpeggioContactProvider> engine;
        result.contacts = (contact_type_ == InteractionType::None)
          ? engine.compute(*result.topology)
          : engine.compute(*result.topology, contact_type_);
      } else {
        InteractionEngine<MolStarContactProvider> engine;
        result.contacts = (contact_type_ == InteractionType::None)
          ? engine.compute(*result.topology)
          : engine.compute(*result.topology, contact_type_);
      }

      result.num_contacts = result.contacts.size();
      result.success = true;

    } catch (const std::exception& e) {
      result.success = false;
      Logger::get_logger()->error("Error processing file {} for contacts: {}", file_path, e.what());
    }

    return result;
  }

private:
  ContactProvider provider_;
  InteractionType contact_type_;
};

} // namespace lahuta::tasks

#endif // LAHUTA_PIPELINE_TASKS_CONTACTS_TASK_HPP
