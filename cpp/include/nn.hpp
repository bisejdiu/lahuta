#ifndef LAHUTA_NN_HPP
#define LAHUTA_NN_HPP

#include "contacts/features.hpp"
#include "nsgrid.hpp"
#include "rings.hpp"
#include <cstdint>
#include <vector>

namespace lahuta {

class Luni;

enum class InteractionType {
  None,
  Any,
  Hydrophobic,
  Halogen,
  HydrogenBond,
  WeakHydrogenBond,
  Ionic,
  MetalCoordination,
  CationPi,
  PiStackingP,
  PiStackingT,
};

enum class EntityType : uint8_t {
  Atom = 0,
  Ring = 1,
  Group = 2,
};

inline std::string entity_type_to_string(EntityType type) {
  switch (type) {
  case EntityType::Atom:
    return "Atom";
  case EntityType::Ring:
    return "Ring";
  case EntityType::Group:
    return "Group";
  default:
    return "Unknown";
  }
}

// Number of bits to shift the entity type into the upper 8 bits of a 64-bit value.
constexpr uint64_t ENTITY_TYPE_SHIFT = 56;

// Mask to isolate the lower 56 bits for the entity index.
constexpr uint64_t ENTITY_INDEX_MASK = 0x00FFFFFFFFFFFFFF;

// EntityID packs an EntityType (upper 8 bits) and EntityIndex (lower 56 bits) into a uint64_t.
// - ENTITY_TYPE_SHIFT shifts the EntityType into the uppermost byte.
// - ENTITY_INDEX_MASK ensures the EntityIndex fits within the lower 56 bits.
using EntityID = uint64_t;
using EntityIndex = uint64_t;

//! Create and extract EntityID components
inline EntityID make_entity_id(EntityType type, EntityIndex index) {
  return (static_cast<uint64_t>(type) << ENTITY_TYPE_SHIFT) | (index & ENTITY_INDEX_MASK);
}

//! Extract the EntityType from an EntityID
inline EntityType get_entity_type(EntityID id) { return static_cast<EntityType>(id >> ENTITY_TYPE_SHIFT); }

//! Extract the EntityIndex from an EntityID
inline uint64_t get_entity_index(EntityID id) { return id & ENTITY_INDEX_MASK; }

//! Compare two EntityTypes by priority
inline bool is_higher_priority(EntityType type1, EntityType type2) {
  return static_cast<uint8_t>(type1) < static_cast<uint8_t>(type2);
}

struct EntityVisitor {
  void visit(const RDKit::Atom &atom) const {}
  void visit(const RingData &ring) const {}
  void visit(const Feature &feature) const {}
};

struct Contact {
  EntityID entity1;
  EntityID entity2;
  float distance;
  InteractionType type;

  Contact(EntityID e1, EntityID e2, float d) : Contact(e1, e2, d, InteractionType::Any) {}
  Contact(EntityID e1, EntityID e2, float d, InteractionType t)
      : entity1(e1), entity2(e2), distance(d), type(t) {

    EntityType type1 = get_entity_type(entity1);
    EntityType type2 = get_entity_type(entity2);

    // Sort order:      Atom < Ring < Group
    // different types: sort by type
    // same types:      sort by index
    if (type1 != type2) {
      if (!is_higher_priority(type1, type2)) {
        std::swap(entity1, entity2);
      }
    } else {
      if (entity1 > entity2) {
        std::swap(entity1, entity2);
      }
    }
  }

  bool operator==(const Contact &other) const {
    return (entity1 == other.entity1 && entity2 == other.entity2) ? (type == other.type) : false;
  }

  bool operator<(const Contact &other) const {
    return (entity1 != other.entity1)   ? (entity1 < other.entity1)
           : (entity2 != other.entity2) ? (entity2 < other.entity2)
                                        : (type < other.type);
  }

  bool operator>(const Contact &other) const {
    if (*this == other) {
      return false;
    }
    return !(*this < other);
  }

  bool operator!=(const Contact &other) const { return !(*this == other); }

  bool operator<=(const Contact &other) const {
    if (*this < other || *this == other) {
      return true;
    }
    return false;
  }

  bool operator>=(const Contact &other) const {
    if (*this > other || *this == other) {
      return true;
    }
    return false;
  }
};

// Container for storing interactions (Interaction Container - IC)
class Contacts {
public:
  std::vector<Contact> interactions;
  bool is_sorted = false;
  std::string instance_name{};

  Contacts() = default;
  Contacts(const Luni *luni) : luni(luni), is_sorted(false) {}

  void set_luni(const Luni *luni) { this->luni = luni; }

  void sort_interactions() {
    std::sort(interactions.begin(), interactions.end());
    is_sorted = true;
  }

  void sort_if_not_sorted() {
    if (!is_sorted) {
      sort_interactions();
    }
  }

  static void prepare_input(Contacts &lhs, Contacts &rhs) {
    lhs.sort_if_not_sorted();
    rhs.sort_if_not_sorted();
  }

  void add(Contacts &ic) { add(ic.interactions); }

  void add(const Contact &interaction) {
    interactions.push_back(interaction);
    is_sorted = false;
  }

  void add(EntityID e1, EntityID e2, float d, InteractionType t) { add(Contact(e1, e2, d, t)); }

  void add(std::vector<Contact> &interactions) {
    this->interactions.insert(this->interactions.end(), interactions.begin(), interactions.end());
    is_sorted = false;
  }

  void add_many(
      const NSResults &neighbors, const std::vector<EntityID> &e1, const std::vector<EntityID> &e2,
      InteractionType type = InteractionType::Any);

  void add_many(
      const NSResults &neighbors, const std::vector<EntityID> &entities,
      InteractionType type = InteractionType::Any) {
    add_many(neighbors, entities, entities, type);
  }

  // FIX: The visit_entity functions represent a proof of concept for calling 
  // functions and passing the correct type of object to the function.
  void visit_entity(const Luni &luni, EntityID entity) const;

  template <typename Func1, typename Func2>
  void visit_entity(const Luni &luni, EntityID entity, Func1 func1, Func2 func2) const;

  template <typename Func> void visit_entity(const Luni &luni, EntityID entity, Func func) const {
    visit_entity(luni, entity, func, func);
  }

  template <typename Func> void visit_entities(Func func) const {
    for (const auto &interaction : interactions) {
      visit_entity(*luni, interaction.entity1, func);
      visit_entity(*luni, interaction.entity2, func);
    }
  }

  // Make all interactions generic (i.e. remove type information)
  // Leads to loss of information.
  void make_generic();

  /// union of two Contacts
  template <typename T> Contacts &set_union(T &&other);

  /// intersection of two Contacts
  template <typename T> Contacts set_intersection(T &&other);

  /// difference of two Contacts
  template <typename T> Contacts set_difference(T &&other);

  /// symmetric difference of two Contacts
  template <typename T> Contacts set_symmetric_difference(T &&other);

  bool operator==(Contacts &other);

  bool operator!=(Contacts &other) { return !(*this == other); }

  /// intersection of two Contacts: C1 & C2
  Contacts operator&(Contacts &other) { return set_intersection(other); }

  /// union of two Contacts: C1 | C2
  Contacts operator|(Contacts &other) { return set_union(other); }

  /// difference of two Contacts: C1 - C2
  Contacts operator-(Contacts &other);

  /// symmetric difference of two Contacts: C1 ^ C2
  Contacts operator^(Contacts &other);

  /// union of C1 with C2: C1 |= C2
  template <typename T> Contacts &operator|=(T &&other);

  /// add C2 to C1: C1 += C2
  template <typename T> Contacts &operator+=(T &&other);

  /// subtract C2 from C1: C1 -= C2
  template <typename T> Contacts &operator-=(T &&other);

  /// intersect C1 with C2: C1 &= C2
  template <typename T> Contacts &operator&=(T &&other);

  /// symmetric difference of C1 with C2: C1 ^= C2
  template <typename T> Contacts &operator^=(T &&other);

  /// access interaction by index
  /// returns a Contacts object with a single interaction
  Contacts operator[](int index) const { return Contacts(luni, interactions[index]); }

  int size() const { return interactions.size(); }
  const Luni *get_luni() const { return luni; }
  void print_interactions() const;

  friend class Luni;

private:
  Contacts(const Luni *luni, std::vector<Contact> interactions) : luni(luni), interactions(interactions) {}
  Contacts(const Luni *luni, Contact interaction) : luni(luni) { interactions.push_back(interaction); }

  std::string get_entity_atoms(const EntityID &entity) const;

private:
  const Luni *luni;
  const EntityVisitor visitor;
};

// clang-format off

inline void Contacts::make_generic() {
  for (auto &interaction : interactions) {
    interaction.type = InteractionType::Any;
  }
  // remove duplicates and re-sort
  std::sort(interactions.begin(), interactions.end());
  interactions.erase(
    std::unique(interactions.begin(),
                interactions.end()),
    interactions.end()
  );
}

template <typename T> Contacts &Contacts::set_union(T &&other) {
  prepare_input(*this, other);
  std::vector<Contact> result;
  result.reserve(interactions.size() + other.interactions.size());
  std::set_union(interactions.begin(),
                 interactions.end(),
                 other.interactions.begin(),
                 other.interactions.end(),
                 std::back_inserter(result)
  );
  interactions = std::move(result);
  is_sorted = true;
  return *this;
}

template <typename T> Contacts Contacts::set_intersection(T &&other) {
  prepare_input(*this, other);
  Contacts result(this->luni);
  std::set_intersection(interactions.begin(),
                        interactions.end(),
                        other.interactions.begin(),
                        other.interactions.end(),
                        std::back_inserter(result.interactions)
  );
  return result;
}

template <typename T> Contacts Contacts::set_difference(T &&other) {
  prepare_input(*this, other);
  Contacts result(this->luni);
  std::set_difference(interactions.begin(),
                      interactions.end(),
                      other.interactions.begin(),
                      other.interactions.end(),
                      std::back_inserter(result.interactions)
  );
  return result;
}

template <typename T> Contacts Contacts::set_symmetric_difference(T &&other) {
  prepare_input(*this, other);

  Contacts result(this->luni);
  std::set_symmetric_difference(interactions.begin(),
                                interactions.end(),
                                other.interactions.begin(),
                                other.interactions.end(),
                                std::back_inserter(result.interactions)
  );

  return result;
}

// clang-format on

inline bool Contacts::operator==(Contacts &other) {
  prepare_input(*this, other);
  return interactions == other.interactions;
}

inline Contacts Contacts::operator-(Contacts &other) { return set_difference(other); }

inline Contacts Contacts::operator^(Contacts &other) {
  Contacts result = set_union(other);
  Contacts common = set_intersection(other);
  return result - common;
}

template <typename T> Contacts &Contacts::operator|=(T &&other) {
  Contacts result = set_union(std::forward<T>(other));
  interactions = std::move(result.interactions);
  is_sorted = true;
  return *this;
}

template <typename T> Contacts &Contacts::operator+=(T &&other) {
  interactions.insert(interactions.end(), other.interactions.begin(), other.interactions.end());
  is_sorted = false;
  return *this;
}

template <typename T> Contacts &Contacts::operator-=(T &&other) {
  Contacts result = set_difference(std::forward<T>(other));
  interactions = std::move(result.interactions);
  is_sorted = true;
  return *this;
}

template <typename T> Contacts &Contacts::operator&=(T &&other) {
  Contacts result = set_intersection(std::forward<T>(other));
  interactions = std::move(result.interactions);
  is_sorted = true;
  return *this;
}

template <typename T> Contacts &Contacts::operator^=(T &&other) {
  Contacts result = set_symmetric_difference(std::forward<T>(other));
  interactions = std::move(result.interactions);
  is_sorted = true;
  return *this;
}

} // namespace lahuta

#endif // LAHUTA_NN_HPP
