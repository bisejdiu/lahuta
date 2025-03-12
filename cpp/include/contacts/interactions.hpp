#ifndef LAHUTA_INTERACTIONS_HPP
#define LAHUTA_INTERACTIONS_HPP

#include "neighbors.hpp"
#include <optional>

namespace lahuta {

class Luni;

class InteractionBase {
  Contacts contacts;

public:
  virtual void handle(const GroupEntity &feature1, const GroupEntity &feature2) = 0;
  virtual ~InteractionBase() = default;
  void add_interaction(const Contact &interaction) { contacts.add(interaction); }
  Contacts release() { return std::move(contacts); }
};

// TODO: 1. We need to provide support for controling contact options
//       2. InteractionOptions is used by the contact functions making compilation slow

struct InteractionOptions {
  float distance_cutoff;
  bool hbond = true;
  bool weak_hbond = true;
  bool hydrophobic = true;
  bool halogen = true;
  bool ionic = true;
  bool metalic = true;
  bool cationpi = true;
  bool pistacking = true;

  bool sort_globally = true;
  bool sort_per_contact_type = false;
};

class Interactions {

  struct InteractionEntry {
    Contacts (Interactions::*func)() const;
    bool InteractionOptions::*flag;
  };

  using InteractionFunc = Contacts (Interactions::*)() const;
  const InteractionEntry entries[8] = {
      {&Interactions::hbond, &InteractionOptions::hbond},
      {&Interactions::weak_hbond, &InteractionOptions::weak_hbond},
      {&Interactions::halogen, &InteractionOptions::halogen},
      {&Interactions::hydrophobic, &InteractionOptions::hydrophobic},
      {&Interactions::metalic, &InteractionOptions::metalic},
      {&Interactions::ionic, &InteractionOptions::ionic},
      {&Interactions::cationpi, &InteractionOptions::cationpi},
      {&Interactions::pistacking, &InteractionOptions::pistacking}};

public:
  Interactions(const Luni &luni, std::optional<InteractionOptions> opts = std::nullopt)
      : luni_(luni), opts_(opts.value_or(InteractionOptions{5.0})) {}

  [[nodiscard]] Contacts hbond() const;
  [[nodiscard]] Contacts weak_hbond() const;
  [[nodiscard]] Contacts hydrophobic() const;
  [[nodiscard]] Contacts halogen() const;
  [[nodiscard]] Contacts ionic() const;
  [[nodiscard]] Contacts metalic() const;
  [[nodiscard]] Contacts cationpi() const;
  [[nodiscard]] Contacts pistacking() const;

  Contacts compute_contacts() {
    Contacts result{&luni_};

    for (const auto &entry : entries) {
      if (!(opts_.*(entry.flag))) continue;
      /*result.add((this->*(entry.func))());*/
      Contacts cts = (this->*(entry.func))();
      if (opts_.sort_per_contact_type) cts.sort_interactions();
      result.add(cts);
    }

    if (opts_.sort_globally) result.sort_interactions();
    std::cout << "Total Contacts: " << result.size() << std::endl;
    return result;
  }

private:
  const Luni &luni_;
  /*const std::vector<Feature> *group_features_;*/
  InteractionOptions opts_{5.0};
};

} // namespace lahuta

#endif // LAHUTA_INTERACTIONS_HPP
