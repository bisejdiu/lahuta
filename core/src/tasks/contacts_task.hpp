#ifndef LAHUTA_PIPELINE_TASKS_CONTACTS_TASK_HPP
#define LAHUTA_PIPELINE_TASKS_CONTACTS_TASK_HPP

#include "contacts/arpeggio/provider.hpp"
#include "contacts/engine.hpp"
#include "contacts/molstar/provider.hpp"
#include "entities/contact.hpp"
#include "entities/interaction_types.hpp"
#include "lahuta.hpp"
#include "logging.hpp"
#include "models/topology.hpp"
#include "tasks/model_db_writer.hpp"
#include "topology.hpp"
#include <chemistry/atom_typing.hpp>
#include <chemistry/group_typing.hpp>
#include <memory>
#include <string>
#include <string_view>
#include <typing/flags.hpp>
#include <valence_model.hpp>

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
    std::unique_ptr<Topology> topology;
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

  // Overload for database-sourced model data
  result_type operator()(const tasks::ModelWriteTask::result_type& model_data) const {
    result_type result;
    result.file_path    = model_data.file_path;
    result.provider     = provider_;
    result.contact_type = contact_type_;
    result.success      = false;
    result.num_contacts = 0;
    result.topology     = nullptr;

    try {

      auto mol = std::make_shared<RDKit::RWMol>();
      lahuta::build_model_topology(mol, model_data.data, ModelTopologyMethod::CSR);

      auto luni = Luni::create(mol);

      const auto computer_type = (provider_ == ContactProvider::Arpeggio)
        ? ContactComputerType::Arpeggio
        : ContactComputerType::Molstar;

      luni.set_atom_typing_method(computer_type);

      result.topology = luni.release_topology();

      auto &residues = result.topology->get_residues();
      residues.build();

      auto &data = result.topology->get_engine().get_data();

      data.atoms.reserve(data.mol->getNumAtoms());
      for (const auto &atom : data.mol->atoms()) {
        auto atom_type = atom->getCompAtomType();
        data.atoms.push_back(AtomRec{static_cast<AtomType>(atom_type), *atom});
      }

      // data.mol->updatePropertyCache(false);
      // ValenceModel valence_model;
      // valence_model.apply(*data.mol);
      // data.atoms  = AtomTypeAnalysis()(*data.mol);

      data.groups = GroupTypeAnalysis::analyze(*data.mol, *data.residues);
      data.rings  = topology::AtomTypingKernel::populate_ring_entities(*data.mol);

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
      Logger::get_logger()->error("Error processing model data for {}: {}", model_data.file_path, e.what());
    }

    return result;
  }

private:
  ContactProvider provider_;
  InteractionType contact_type_;
};

} // namespace lahuta::tasks

#endif // LAHUTA_PIPELINE_TASKS_CONTACTS_TASK_HPP
