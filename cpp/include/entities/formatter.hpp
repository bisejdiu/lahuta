#ifndef LAHUTA_ENTITIES_FORMATTER_HPP
#define LAHUTA_ENTITIES_FORMATTER_HPP

#include "contact.hpp"
#include "entity_id.hpp"
#include "interaction_types.hpp"
#include "records.hpp"
#include "topology.hpp"
#include <GraphMol/MonomerInfo.h>
#include <iostream>
#include <vector>

// clang-format off
namespace lahuta {

class ContactTableFormatter {
private:
  static constexpr size_t MAX_ATOM_NAME_LEN = 3;
  static constexpr size_t MAX_RES_NAME_LEN = 4;
  static constexpr size_t MAX_CHAIN_ID_LEN = 2;
  static constexpr size_t MAX_RING_ATOMS = 8;
  static constexpr size_t FALLBACK_RING_ATOMS = 4;

  static constexpr size_t SEPARATOR_WIDTH = 7; // "  -&-  "
  static constexpr size_t DISTANCE_WIDTH = 8;  // "12.34567"
  static constexpr size_t TYPE_WIDTH = 20;     // interaction type names

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

    // Max width needed for each column
    size_t max_entity1_width = 0;
    size_t max_entity2_width = 0;
    std::vector<std::pair<std::string, std::string>> formatted_entities;
    formatted_entities.reserve(contact_set.size());

    for (const auto& contact : contact_set) {
      std::string entity1_info = format_entity_compact(topology, contact.lhs);
      std::string entity2_info = format_entity_compact(topology, contact.rhs);

      max_entity1_width = std::max(max_entity1_width, entity1_info.length());
      max_entity2_width = std::max(max_entity2_width, entity2_info.length());
      formatted_entities.emplace_back(std::move(entity1_info), std::move(entity2_info));
    }

    // Print with calculated widths for each column
    size_t contact_idx = 0;
    for (const auto& contact : contact_set) {
      const auto& [entity1_info, entity2_info] = formatted_entities[contact_idx++];
      std::string type_str = interaction_type_to_string(contact.type);

      std::cout << std::left
                << std::setw(max_entity1_width) << entity1_info
                << "  -&-  "
                << std::setw(max_entity2_width) << entity2_info
                << " "
                << std::setw(DISTANCE_WIDTH) << std::fixed << std::setprecision(3) << contact.distance
                << " "
                << type_str << std::endl;
    }

    if (!interaction_name.empty()) {
      std::cout << "Total " << interaction_name << " contacts: " << contact_set.size() << std::endl;
      std::cout << std::endl;
    }
  }
};

} // namespace lahuta

#endif // LAHUTA_ENTITIES_FORMATTER_HPP 
