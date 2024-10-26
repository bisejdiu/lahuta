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
    const RingData &ring = _luni.get_entity<RingData>(entity);
    visitor.visit(ring);
    break;
  }
  case EntityType::Group: {
    const Feature &group = _luni.get_entity<Feature>(entity);
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
    const RingData &ring = luni.get_entity<RingData>(entity);
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

    visit_entity(*luni, interaction.entity1);
    visit_entity(*luni, interaction.entity2);

    auto type = get_entity_type(interaction.entity1);
    std::string e1_atoms = "";
    std::string e2_atoms = "";
    if (type == EntityType::Group) {
      const Feature &group = luni->get_entity<Feature>(interaction.entity1);
      for (const auto *atom : group.members) {
        e1_atoms += std::to_string(atom->getIdx()) + " ";
      }
      const Feature &group2 = luni->get_entity<Feature>(interaction.entity2);
      for (const auto *atom : group2.members) {
        e2_atoms += std::to_string(atom->getIdx()) + " ";
      }
    }
    if (type == EntityType::Atom) {
      const RDKit::Atom &atom = luni->get_entity<RDKit::Atom>(interaction.entity1);
      e1_atoms = std::to_string(atom.getIdx());
      const RDKit::Atom &atom2 = luni->get_entity<RDKit::Atom>(interaction.entity2);
      e2_atoms = std::to_string(atom2.getIdx());
    }
    if (type == EntityType::Ring) {
      const RingData &ring = luni->get_entity<RingData>(interaction.entity1);
      auto atom_ids = ring.atom_ids();
      for (const auto &atom_id : atom_ids) {
        e1_atoms += std::to_string(atom_id) + " ";
      }
      const RingData &ring2 = luni->get_entity<RingData>(interaction.entity2);
      auto atom_ids2 = ring2.atom_ids();
      for (const auto &atom_id : atom_ids2) {
        e2_atoms += std::to_string(atom_id) + " ";
      }
    }
    std::cout << "Interaction between entity " << entity_type_to_string(get_entity_type(interaction.entity1))
              << " (" << get_entity_index(interaction.entity1) << ") and entity "
              << entity_type_to_string(get_entity_type(interaction.entity2)) << " ("
              << get_entity_index(interaction.entity2) << ") with distance " << interaction.distance
              << " --- " << e1_atoms << " --- " << e2_atoms << "\n";
  }
}

} // namespace lahuta
