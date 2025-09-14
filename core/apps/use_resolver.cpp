#include <contacts/engine.hpp>
#include <contacts/molstar/provider.hpp>
#include <iostream>
#include <sstream>
#include <string>
#include <variant>

#include "entities/interaction_types.hpp"
#include "entities/resolver.hpp"
#include "lahuta.hpp"
#include <GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RDKitBase.h>

using namespace lahuta;

// clang-format off
static std::string atom_tuple(const RDKit::Atom &atom) {

  const auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom.getMonomerInfo());
  if (!info) return std::to_string(atom.getIdx()) + "-UNK-UNK-0";

  const int res_num = info->getResidueNumber();
  const std::string atom_name = info->getName();
  const std::string res_name  = info->getResidueName();
  std::ostringstream oss;
  oss << atom.getIdx() << '-' << atom_name << '-' << res_name << '-' << res_num;
  return oss.str();
}

static std::string describe_entity(const EntityRef &v) {
  return std::visit(
      [](const auto &ref) -> std::string {
        using RefT = std::decay_t<decltype(ref)>;
        using RawT = std::remove_cv_t<typename RefT::type>; // reference_wrapper<T>, strip const
        const auto &rec = ref.get();
        if constexpr (std::is_same_v<RawT, AtomRec>) {
          return atom_tuple(rec.atom.get());
        } else if constexpr (std::is_same_v<RawT, RingRec>) {
          std::ostringstream oss;
          oss << '(';
          for (size_t i = 0; i < rec.atoms.size(); ++i) {
            if (i > 0) oss << ", ";
            oss << atom_tuple(rec.atoms[i].get());
          }
          oss << ')';
          return oss.str();
        } else {
          static_assert(std::is_same_v<RawT, GroupRec>, "Unhandled record type");
          std::ostringstream oss;
          oss << '(';
          for (size_t i = 0; i < rec.atoms.size(); ++i) {
            if (i > 0) oss << ", ";
            oss << atom_tuple(rec.atoms[i].get());
          }
          oss << ')';
          return oss.str();
        }
      },
      v);
}

int main(int argc, char **argv) {
  if (argc != 2 || std::string(argv[1]) == "--help") {
    std::cout << "Usage: use_resolver input_file.cif\n";
    return argc == 2 ? 0 : 1;
  }

  Luni luni(argv[1]);
  if (!luni.build_topology()) {
    std::cerr << "Failed to build topology from " << argv[1] << "\n";
    return 1;
  }
  const Topology &top = luni.get_topology();

  InteractionEngine<MolStarContactProvider> engine;
  ContactSet contacts = engine.compute(top);

  EntityResolver resolver(top);
  const auto resolved = resolver.resolve_all(contacts);

  const std::size_t limit = std::min<std::size_t>(resolved.size(), 10000);
  for (std::size_t i = 0; i < limit; ++i) {
    const auto &c = contacts.data()[i];
    const auto &pair = resolved[i];
    std::cout << describe_entity(pair.first) << " " << describe_entity(pair.second) << " "
              << interaction_type_to_string(c.type) << " " << c.distance << "\n";
  }

  return 0;
}
