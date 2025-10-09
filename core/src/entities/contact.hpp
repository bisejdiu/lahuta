#ifndef LAHUTA_ENTITIES_CONTACT_HPP
#define LAHUTA_ENTITIES_CONTACT_HPP

#include <vector>
#include <algorithm>

#include "entity_id.hpp"
#include "interaction_types.hpp"

// clang-format off
namespace lahuta {

//
// Contact represents an interaction between two entities.
// Always stored with lhs < rhs for consistent ordering.
//
struct Contact {
  EntityID lhs;
  EntityID rhs;
  float distance;
  InteractionType type;

  Contact(EntityID e1, EntityID e2, float dist, InteractionType t);

  bool operator==(const Contact& other) const {
    return lhs.raw == other.lhs.raw &&  rhs.raw == other.rhs.raw &&  type == other.type;
  }

  bool operator!=(const Contact& other) const {
    return !(*this == other);
  }

  bool operator<(const Contact& other) const {
    return std::tie(lhs.raw, rhs.raw, type) <  std::tie(other.lhs.raw, other.rhs.raw, other.type);
  }
};

struct ContactHash { // FIX: used?
  std::size_t operator()(const Contact& c) const {
    std::size_t h1 = std::hash<uint64_t>{}(c.lhs.raw);
    std::size_t h2 = std::hash<uint64_t>{}(c.rhs.raw);
    std::size_t h3 = std::hash<uint32_t>{}(static_cast<uint32_t>(c.type));
    return (h1 ^ (h2 << 1)) ^ (h3 << 2);
  }
};

//
// ContactSet is an always-ordered container for Contact objects.
// It provides efficient set operations and guarantees no duplicates.
//
class ContactSet {
public:
  ContactSet() = default;
  ContactSet(std::initializer_list<Contact> contacts);
  template<typename Vt, typename = std::enable_if_t<std::is_same<std::decay_t<Vt>, std::vector<Contact>>::value>>
  explicit ContactSet(Vt&& v, bool make_unique = false) : contacts_(std::forward<Vt>(v)) {
    ensure_sorted();

    if (!make_unique) return;
    contacts_.erase(std::unique(contacts_.begin(), contacts_.end()), contacts_.end());
  }

  void insert(const Contact& contact);
  void insert(const ContactSet& other);

  ContactSet set_union(const ContactSet& other) const;
  ContactSet set_intersection(const ContactSet& other) const;
  ContactSet set_difference(const ContactSet& other) const;
  ContactSet set_symmetric_difference(const ContactSet& other) const;

  ContactSet operator|(const ContactSet& other) const { return set_union(other); }
  ContactSet operator&(const ContactSet& other) const { return set_intersection(other); }
  ContactSet operator-(const ContactSet& other) const { return set_difference(other); }
  ContactSet operator^(const ContactSet& other) const { return set_symmetric_difference(other); }

  ContactSet& operator|=(const ContactSet& other);
  ContactSet& operator&=(const ContactSet& other);
  ContactSet& operator-=(const ContactSet& other);
  ContactSet& operator^=(const ContactSet& other);

  void make_generic();

  auto begin() const { return contacts_.begin(); }
  auto end() const { return contacts_.end(); }

  bool empty() const { return contacts_.empty(); }
  size_t size() const { return contacts_.size(); }
  const auto& data() const { return contacts_; }

  const Contact& operator[](size_t index) const { return contacts_[index]; }
  Contact& operator[](size_t index) { return contacts_[index]; }

private:
    std::vector<Contact> contacts_;

    void ensure_sorted();
};

} // namespace lahuta

#endif // LAHUTA_ENTITIES_CONTACT_HPP 
