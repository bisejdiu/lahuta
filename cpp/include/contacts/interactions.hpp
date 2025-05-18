#ifndef LAHUTA_INTERACTIONS_HPP
#define LAHUTA_INTERACTIONS_HPP

#include "entities/contact.hpp"

namespace lahuta {

class Topology;

struct InteractionOptions {
  float distance_cutoff = 5.0;

  bool hbond       = true;
  bool weak_hbond  = true;
  bool hydrophobic = true;
  bool halogen     = true;
  bool ionic       = true;
  bool metal       = true;
  bool cationpi    = true;
  bool pistacking  = true;

  bool sort_globally = true;
};

class Interactions {
public:

  Interactions(const Topology &topology, InteractionOptions opts = InteractionOptions{})
      : topology_(topology), opts_(opts) {}

  [[nodiscard]] ContactSet hbond() const;
  [[nodiscard]] ContactSet weak_hbond() const;
  [[nodiscard]] ContactSet hydrophobic() const;
  [[nodiscard]] ContactSet halogen() const;
  [[nodiscard]] ContactSet ionic() const;
  [[nodiscard]] ContactSet metal() const;
  [[nodiscard]] ContactSet cationpi() const;
  [[nodiscard]] ContactSet pistacking() const;

  ContactSet compute_contacts() {
    ContactSet result;

    if (opts_.hbond) {
      auto contacts = hbond();
      result |= contacts;
    }

    if (opts_.weak_hbond) {
      auto contacts = weak_hbond();
      result |= contacts;
    }

    if (opts_.hydrophobic) {
      auto contacts = hydrophobic();
      result |= contacts;
    }

    if (opts_.halogen) {
      auto contacts = halogen();
      result |= contacts;
    }

    if (opts_.ionic) {
      auto contacts = ionic();
      result |= contacts;
    }

    if (opts_.metal) {
      auto contacts = metal();
      result |= contacts;
    }

    if (opts_.cationpi) {
      auto contacts = cationpi();
      result |= contacts;
    }

    if (opts_.pistacking) {
      auto contacts = pistacking();
      result |= contacts;
    }

    return result;
  }

private:
  const Topology &topology_;
  InteractionOptions opts_;
};

} // namespace lahuta

#endif // LAHUTA_INTERACTIONS_HPP
