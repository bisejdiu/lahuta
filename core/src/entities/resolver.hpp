#ifndef LAHUTA_ENTITIES_RESOLVER_HPP
#define LAHUTA_ENTITIES_RESOLVER_HPP

#include <functional>
#include <utility>
#include <variant>
#include <vector>

#include "contact.hpp"
#include "entity_id.hpp"
#include "records.hpp"
#include "topology.hpp"

// clang-format off
namespace lahuta {

using EntityRef = std::variant<
  std::reference_wrapper<const AtomRec>,
  std::reference_wrapper<const RingRec>,
  std::reference_wrapper<const GroupRec>
>;

struct EntityResolver {
  explicit EntityResolver(const Topology &top) : top_(&top) {}

  const Topology &topology() const noexcept { return *top_; }

  EntityRef resolve(EntityID eid) const {
    switch (eid.kind()) {
      case Kind::Atom:  return std::cref(top_->resolve<Kind::Atom> (eid));
      case Kind::Ring:  return std::cref(top_->resolve<Kind::Ring> (eid));
      case Kind::Group: return std::cref(top_->resolve<Kind::Group>(eid));
      default:          throw std::invalid_argument("Invalid EntityID kind");
    }
  }

  std::pair<EntityRef, EntityRef> resolve(const Contact &c) const {
    return {resolve(c.lhs), resolve(c.rhs)};
  }

  std::vector<std::pair<EntityRef, EntityRef>> resolve_all(const ContactSet &set) const {
    // no deep copies, but it does materialize a container of references.
    std::vector<std::pair<EntityRef, EntityRef>> out;
    out.reserve(set.size());
    for (const auto &c : set.data()) {
      out.emplace_back(resolve(c));
    }
    return out;
  }

  // Visit a single entity with a callable that accepts any of AtomRec, RingRec, GroupRec
  template <class Fn>
  decltype(auto) visit(EntityID id, Fn&& fn) const {
    switch (id.kind()) {
      case Kind::Atom:  return std::forward<Fn>(fn)(top_->resolve<Kind::Atom>(id));
      case Kind::Ring:  return std::forward<Fn>(fn)(top_->resolve<Kind::Ring>(id));
      case Kind::Group: return std::forward<Fn>(fn)(top_->resolve<Kind::Group>(id));
      default: throw std::runtime_error("Invalid EntityID kind");
    }
  }

  // Visit both sides of a contact
  template <class Fn>
  void visit(const Contact& c, Fn&& fn) const {
    visit(c.lhs, [&](const auto& a){ visit(c.rhs, [&](const auto& b){ fn(a,b); }); });
  }

  // Iterate contacts without allocating an intermediate container
  template <class Fn>
  void for_each(const ContactSet& set, Fn&& fn) const {
    for (const auto& c : set) visit(c, fn);
  }

private:
  const Topology *top_;
};

} // namespace lahuta

#endif // LAHUTA_ENTITIES_RESOLVER_HPP
