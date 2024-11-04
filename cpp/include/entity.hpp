#ifndef LAHUTA_ENTITY_HPP
#define LAHUTA_ENTITY_HPP

#include "GraphMol/Atom.h"
#include <cstdint>
#include <string>

namespace lahuta {

enum class EntityType : uint8_t {
  Atom = 0,
  Ring = 1,
  Group = 2,
};

/// Number of bits to shift the entity type into the upper 8 bits of a 64-bit value.
constexpr uint64_t ENTITY_TYPE_SHIFT = 56;

/// Mask to isolate the lower 56 bits for the entity index.
constexpr uint64_t ENTITY_INDEX_MASK = 0x00FFFFFFFFFFFFFF;

// EntityID packs an EntityType (upper 8 bits) and EntityIndex (lower 56 bits) into a uint64_t.
// - ENTITY_TYPE_SHIFT shifts the EntityType into the uppermost byte.
// - ENTITY_INDEX_MASK ensures the EntityIndex fits within the lower 56 bits.
using EntityID = uint64_t;
using EntityIndex = uint64_t;

/// Create and extract EntityID components
inline EntityID make_entity_id(EntityType type, EntityIndex index) {
  return (static_cast<uint64_t>(type) << ENTITY_TYPE_SHIFT) | (index & ENTITY_INDEX_MASK);
}

/// Extract the EntityType from an EntityID
inline EntityType get_entity_type(EntityID id) { return static_cast<EntityType>(id >> ENTITY_TYPE_SHIFT); }

/// Extract the EntityIndex from an EntityID
inline uint64_t get_entity_index(EntityID id) { return id & ENTITY_INDEX_MASK; }

/// Compare two EntityTypes by priority
inline bool is_higher_priority(EntityType type1, EntityType type2) {
  return static_cast<uint8_t>(type1) < static_cast<uint8_t>(type2);
}

class RingEntity;
class GroupEntity;

struct EntityVisitor {
  void visit(const RDKit::Atom &atom) const {}
  void visit(const RingEntity &ring) const {}
  void visit(const GroupEntity &feature) const {}
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

} // namespace lahuta

#endif // LAHUTA_ENTITY_HPP
