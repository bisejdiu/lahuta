#include <algorithm>

#include "entities/contact.hpp"

namespace lahuta {

void ContactSet::ensure_sorted() { std::sort(contacts_.begin(), contacts_.end()); }

Contact::Contact(EntityID e1, EntityID e2, float dist, InteractionType t) : distance(dist), type(t) {
  // lhs < rhs is the canon ordering
  if (e1.raw < e2.raw) {
    lhs = e1;
    rhs = e2;
  } else {
    lhs = e2;
    rhs = e1;
  }
}

ContactSet::ContactSet(std::initializer_list<Contact> contacts) : contacts_(contacts)  { ensure_sorted(); }

void ContactSet::insert(const Contact& contact) {
  contacts_.push_back(contact);
  ensure_sorted();
}

void ContactSet::insert(const ContactSet& other) {
  contacts_.reserve(contacts_.size() + other.contacts_.size());
  contacts_.insert(contacts_.end(), other.contacts_.begin(), other.contacts_.end());

  ensure_sorted();
}

ContactSet ContactSet::set_union(const ContactSet& other) const {
  ContactSet result;
  result.contacts_.reserve(contacts_.size() + other.contacts_.size());

  std::set_union(
          contacts_.begin(),       contacts_.end(),
    other.contacts_.begin(), other.contacts_.end(),
    std::back_inserter(result.contacts_)
  );

  return result;
}

ContactSet ContactSet::set_intersection(const ContactSet& other) const {
  ContactSet result;
  result.contacts_.reserve(std::min(contacts_.size(), other.contacts_.size()));

  std::set_intersection(
          contacts_.begin(),       contacts_.end(),
    other.contacts_.begin(), other.contacts_.end(),
    std::back_inserter(result.contacts_)
  );

  return result;
}

ContactSet ContactSet::set_difference(const ContactSet& other) const {
  ContactSet result;
  result.contacts_.reserve(contacts_.size());

  std::set_difference(
          contacts_.begin(),       contacts_.end(),
    other.contacts_.begin(), other.contacts_.end(),
    std::back_inserter(result.contacts_)
  );

  return result;
}

ContactSet ContactSet::set_symmetric_difference(const ContactSet& other) const {
  ContactSet result;
  result.contacts_.reserve(contacts_.size() + other.contacts_.size());

  std::set_symmetric_difference(
          contacts_.begin(),       contacts_.end(),
    other.contacts_.begin(), other.contacts_.end(),
    std::back_inserter(result.contacts_)
  );

  return result;
}

ContactSet& ContactSet::operator|=(const ContactSet& other) {
  *this = set_union(other);
  return *this;
}

ContactSet& ContactSet::operator&=(const ContactSet& other) {
  *this = set_intersection(other);
  return *this;
}

ContactSet& ContactSet::operator-=(const ContactSet& other) {
  *this = set_difference(other);
  return *this;
}

ContactSet& ContactSet::operator^=(const ContactSet& other) {
  *this = set_symmetric_difference(other);
  return *this;
}

void ContactSet::make_generic() {
  for (auto& contact : contacts_) {
    contact.type = InteractionType::Generic;
  }

  // Re-sort and remove duplicates
  ensure_sorted();
  contacts_.erase(
    std::unique(contacts_.begin(), contacts_.end()),
    contacts_.end()
  );
}

} // namespace lahuta 
