#ifndef LAHUTA_ENTITIES_FORMATTER_HPP
#define LAHUTA_ENTITIES_FORMATTER_HPP

#include "contact.hpp"
#include "entity_id.hpp"
#include "interaction_types.hpp"
#include "records.hpp"
#include "topology.hpp"
#include <GraphMol/MonomerInfo.h>
#include <span.hpp>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

// clang-format off
namespace lahuta {

namespace {
// custom sorting for contacts
static std::vector<Contact> sort_interactions(const ContactSet& contact_set) {
  std::vector<Contact> sorted_contacts(contact_set.begin(), contact_set.end());

  std::stable_sort(sorted_contacts.begin(), sorted_contacts.end(),
    [](const Contact& a, const Contact& b) {
      // compare by type (category only)
      if (a.type.category != b.type.category) {
        return static_cast<uint8_t>(a.type.category) < static_cast<uint8_t>(b.type.category);
      }
      return a < b; // if types are the same, we will use existing ordering
    });

  return sorted_contacts;
}

static std::vector<span<const Contact>> slice_by_type(const std::vector<Contact>& contacts) noexcept {
  std::vector<span<const Contact>> slices;
  const std::size_t n = contacts.size();
  if (n == 0) return slices;

  std::size_t start = 0;
  auto        cur   = contacts[0].type;

  for (std::size_t i = 1; i < n; ++i) {
    if (contacts[i].type != cur) {
      slices.emplace_back(contacts.data() + start, i - start);
      start = i;
      cur   = contacts[i].type;
    }
  }

  slices.emplace_back(contacts.data() + start, n - start); // last slice
  return slices;
}
} // namespace

class ContactTableFormatter {
private:
  static constexpr size_t MAX_ATOM_NAME_LEN   = 3;
  static constexpr size_t MAX_RES_NAME_LEN    = 4;
  static constexpr size_t MAX_CHAIN_ID_LEN    = 2;
  static constexpr size_t MAX_RING_ATOMS      = 8;
  static constexpr size_t FALLBACK_RING_ATOMS = 4;

  static constexpr size_t SEPARATOR_WIDTH = 7;  // "  -&-  "
  static constexpr size_t DISTANCE_WIDTH  = 8;  // "12.34567"
  static constexpr size_t TYPE_WIDTH      = 20; // interaction type names

public:
  static std::string format_truncated_atoms(const std::vector<std::reference_wrapper<const RDKit::Atom>>& atoms) {
    if (atoms.size() <= MAX_RING_ATOMS) {
      std::string result = "(";
      for (size_t i = 0; i < atoms.size(); ++i) {
        if (i > 0) result += ", ";
        const auto& atom = atoms[i].get();
        auto* info = static_cast<const RDKit::AtomPDBResidueInfo*>(atom.getMonomerInfo());
        std::string atom_name = info ? info->getName() : "UNK";
        result += std::to_string(atom.getIdx()) + "-" + atom_name;
      }
      result += ")";
      return result;
    } else {
      // Show first FALLBACK_RING_ATOMS atoms, then "..."
      std::string result = "(";
      for (size_t i = 0; i < FALLBACK_RING_ATOMS; ++i) {
        if (i > 0) result += ", ";
        const auto& atom = atoms[i].get();
        auto* info = static_cast<const RDKit::AtomPDBResidueInfo*>(atom.getMonomerInfo());
        std::string atom_name = info ? info->getName() : "UNK";
        result += std::to_string(atom.getIdx()) + "-" + atom_name;
      }
      result += ", ...)";
      return result;
    }
  }

  static std::string format_entity_compact(const Topology& topology, EntityID entity_id) {
    switch (entity_id.kind()) {
      case Kind::Atom: {
        const auto& atom_rec = topology.atom(entity_id.index());
        const auto& atom = atom_rec.atom.get();
        auto* info = static_cast<const RDKit::AtomPDBResidueInfo*>(atom.getMonomerInfo());
        if (!info) return std::to_string(atom.getIdx()) + "-UNK-0-UNK-X";

        int         res_num   = info->getResidueNumber();
        std::string atom_name = info->getName();
        std::string res_name  = info->getResidueName();
        std::string chain_id  = info->getChainId();

        return std::to_string(atom.getIdx()) + "-" + atom_name + "-" + std::to_string(res_num) + "-" + res_name + "-" + chain_id;
      }
      case Kind::Ring: {
        const auto& ring_rec = topology.ring(entity_id.index());
        if (ring_rec.atoms.empty()) return "()-0-UNK-X";

        std::string atom_info = format_truncated_atoms(ring_rec.atoms);
        const auto& first_atom = ring_rec.atoms[0].get();
        auto* info = static_cast<const RDKit::AtomPDBResidueInfo*>(first_atom.getMonomerInfo());
        if (!info) return atom_info + "-0-UNK-X";

                 int res_num = info->getResidueNumber();
        std::string res_name = info->getResidueName();
        std::string chain_id = info->getChainId();

        return atom_info + "-" + std::to_string(res_num) + "-" + res_name + "-" + chain_id;
      }
      case Kind::Group: {
        const auto& group_rec = topology.group(entity_id.index());
        if (group_rec.atoms.empty()) return "()-0-UNK-X";

        if (group_rec.atoms.size() == 1) {
          // if we only have one atom, we should treat it as an atom entity
          const auto& atom = group_rec.atoms[0].get();
          auto* info = static_cast<const RDKit::AtomPDBResidueInfo*>(atom.getMonomerInfo());
          if (!info) return std::to_string(atom.getIdx()) + "-UNK-0-UNK-X";

                  int res_num   = info->getResidueNumber();
          std::string atom_name = info->getName();
          std::string res_name  = info->getResidueName();
          std::string chain_id  = info->getChainId();

          return std::to_string(atom.getIdx()) + "-" + atom_name + "-" +  std::to_string(res_num) + "-" + res_name + "-" + chain_id;
        }

        std::string atom_info = format_truncated_atoms(group_rec.atoms);

        const auto& first_atom = group_rec.atoms[0].get();
        auto* info = static_cast<const RDKit::AtomPDBResidueInfo*>(first_atom.getMonomerInfo());
        if (!info) return atom_info + "-0-UNK-X";

                 int res_num = info->getResidueNumber();
        std::string res_name = info->getResidueName();
        std::string chain_id = info->getChainId();

        return atom_info + "-" + std::to_string(res_num) + "-" + res_name + "-" + chain_id;
      }
      default:
        return "UNKNOWN-ENTITY";
    }
  }

  static void print_contact_table(const ContactSet& contact_set, const Topology& topology, const std::string& interaction_name = "") {

    auto sorted = sort_interactions(contact_set);
    auto spans  = slice_by_type(sorted);

    for (auto s : spans) {
      InteractionType type = s[0].type;

      // 1) compute widths & buffer formatted entity strings
      std::size_t w1 = 0, w2 = 0;
      std::vector<std::pair<std::string,std::string>> buf;
      buf.reserve(s.size());

      for (auto const& c : s) {
        auto e1 = format_entity_compact(topology, c.lhs);
        auto e2 = format_entity_compact(topology, c.rhs);
        w1 = std::max(w1, e1.size());
        w2 = std::max(w2, e2.size());
        buf.emplace_back(std::move(e1), std::move(e2));
      }

      for (std::size_t i = 0; i < s.size(); ++i) {
        const auto& [e1,e2] = buf[i];
        const auto& c       = s[i];

        std::cout << std::left
                  << std::setw(w1) << e1 << "  -&-  "
                  << std::setw(w2) << e2 << " "
                  << std::setw(DISTANCE_WIDTH) << std::fixed << std::setprecision(3) << c.distance << " "
                  << interaction_type_to_string(type)
                  << "\n";
      }
    }

    if (!interaction_name.empty()) {
      std::cout << "Total " << interaction_name << " contacts: " << contact_set.size() << std::endl;
      std::cout << std::endl;
    }
  }
};

} // namespace lahuta

#endif // LAHUTA_ENTITIES_FORMATTER_HPP 
