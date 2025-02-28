#include "nn.hpp"
#include "lahuta.hpp"

namespace lahuta {

void Contacts::add_many(
    const NSResults &neighbors, const std::vector<EntityID> &e1, const std::vector<EntityID> &e2,
    InteractionType type) {
  for (size_t i = 0; i < neighbors.size(); ++i) {
    auto pair = neighbors.get_pairs()[i];
    auto dist = neighbors.get_distances()[i];
    auto rr = e1[pair.first];
    add(Contact(e1[pair.first], e2[pair.second], dist, type));
  }
  is_sorted = false;
}

void Contacts::visit_entity(const Luni &_luni, EntityID entity) const {
  EntityType type = get_entity_type(entity);

  switch (type) {
  case EntityType::Atom: {
    const RDKit::Atom &atom = _luni.get_entity<RDKit::Atom>(entity);
    visitor.visit(atom);
    break;
  }
  case EntityType::Ring: {
    const RingEntity &ring = _luni.get_entity<RingEntity>(entity);
    visitor.visit(ring);
    break;
  }
  case EntityType::Group: {
    const GroupEntity &group = _luni.get_entity<GroupEntity>(entity);
    visitor.visit(group);
    break;
  }

  default:
    throw std::runtime_error("Unknown entity type");
  }
}

template <typename Func1, typename Func2>
void Contacts::visit_entity(const Luni &luni, EntityID entity, Func1 func1, Func2 func2) const {
  EntityType type = get_entity_type(entity);

  switch (type) {
  case EntityType::Atom: {
    const RDKit::Atom &atom = luni.get_entity<RDKit::Atom>(entity);
    func1(atom);
    break;
  }
  case EntityType::Ring: {
    const RingEntity &ring = luni.get_entity<RingEntity>(entity);
    func2(ring);
    break;
  }
  case EntityType::Group: {
    std::cout << "We will get the feature needed here... " << std::endl;
  }
  default:
    throw std::runtime_error("Unknown entity type");
  }
}

void Contacts::print_interactions() const {
  for (const auto &interaction : interactions) {

    // Get entity information as a formatted string
    std::string entity1_atoms = get_entity_atoms(interaction.entity1);
    std::string entity2_atoms = get_entity_atoms(interaction.entity2);

    // Print formatted interaction details
    std::cout << "Interaction between entity "
              << entity_type_to_string(get_entity_type(interaction.entity1)) << " ("
              << get_entity_index(interaction.entity1) << ") and entity "
              << entity_type_to_string(get_entity_type(interaction.entity2)) << " ("
              << get_entity_index(interaction.entity2) << ") with distance "
              << interaction.distance << " --- " << entity1_atoms << " --- " << entity2_atoms << " "
              << interaction_type_to_string(interaction.type) << "\n";
  }
}

std::string Contacts::get_entity_atoms(const EntityID &entity) const {
  std::string atoms_info;

  switch (get_entity_type(entity)) {
    case EntityType::Group: {
      const GroupEntity &group = luni->get_entity<GroupEntity>(entity);
      for (const auto *atom : group.atoms) {
        atoms_info += std::to_string(atom->getIdx()) + " ";
      }
      break;
    }
    case EntityType::Atom: {
      const RDKit::Atom &atom = luni->get_entity<RDKit::Atom>(entity);
      atoms_info = std::to_string(atom.getIdx());
      break;
    }
    case EntityType::Ring: {
      const RingEntity &ring = luni->get_entity<RingEntity>(entity);
      for (const auto &atom : ring.atoms) {
        atoms_info += std::to_string(atom->getIdx()) + " ";
      }
      break;
    }
    default:
      throw std::runtime_error("Unsupported entity type in get_entity_atoms");
  }

  return atoms_info;
}




} // namespace lahuta
