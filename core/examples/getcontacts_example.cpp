#include <iostream>
#include <string>

#include "compute/topology_snapshot.hpp"
#include "contacts/engine.hpp"
#include "contacts/getcontacts/provider.hpp"
#include "entities/formatter.hpp"
#include "entities/interaction_types.hpp"
#include "entities/resolver.hpp"
#include "lahuta.hpp"

using namespace lahuta;
namespace C = lahuta::compute;

int main(int argc, char **argv) {
  if (argc != 2 || std::string(argv[1]) == "--help") {
    std::cout << "Usage: getcontacts_example <structure.cif>\n";
    return (argc == 2) ? 0 : 1;
  }

  Luni luni(argv[1]);
  if (!luni.build_topology()) {
    std::cerr << "Failed to build topology from " << argv[1] << "\n";
    return 1;
  }

  const auto *topo_ptr = luni.get_topology().get();
  if (!topo_ptr) {
    std::cerr << "Topology is null after build_topology().\n";
    return 1;
  }

  auto &topo = const_cast<Topology &>(*topo_ptr);
  topo.assign_typing(AtomTypingMethod::GetContacts);

  InteractionEngine<GetContactsProvider> engine;
  auto snapshot       = C::snapshot_of(topo, topo.conformer());
  ContactSet contacts = engine.compute(snapshot);

  std::cout << "Computed " << contacts.size() << " contacts via GetContacts provider." << std::endl;

  EntityResolver resolver(topo);
  const auto resolved = resolver.resolve_all(contacts);
  for (std::size_t i = 0; i < resolved.size(); ++i) {
    const auto &pair    = resolved[i];
    const auto &contact = contacts.data()[i];
    std::cout << ContactTableFormatter::format_entity_compact(topo, pair.first) << " \t"
              << ContactTableFormatter::format_entity_compact(topo, pair.second) << " \t"
              << interaction_type_to_string(contact.type) << " \t" << contact.distance << '\n';
  }

  return 0;
}
