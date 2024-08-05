#ifndef LAHUTA_BONDS_HPP
#define LAHUTA_BONDS_HPP

#include <string>

#include "rules.hpp"
#include "token-gperf-generated.hpp"

#include "GraphMol/Atom.h"
#include "GraphMol/MonomerInfo.h"

using namespace lahuta;

// Function pointer type for bond rules
using HandlerFunc = int (*)(std::string_view s2, std::string_view s3);

constexpr inline auto default_fn = [](std::string_view s2, std::string_view s3) { return 1; };

constexpr inline int default_amino_acid_fn(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::default_amino_acid) {
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
  return 1;
}

constexpr inline int default_base_fn(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::default_base) {
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
  return 1;
}

constexpr inline int his_fn(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::his_rules) {
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
    return 1;
}

constexpr inline int arg(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::arg_rule) {
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
    return 1;
}

constexpr inline int phe_fn(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::phe_rule) {
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
    return 1;
}

constexpr inline int trp_fn(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::trp_rule) { 
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
    return 1;
}

constexpr inline int asn(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::asn_rule) { 
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
    return 1;
}

constexpr inline int gln(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::gln_rule) {
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
    return 1;
}

constexpr inline int tyr_fn(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::tyr_rule) { 
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
    return 1;
}

constexpr inline int asp_fn(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::asp_rule) {
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
    return 1;
}

constexpr inline int glu_fn(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::glu_rule) {
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
    return 1;
}

constexpr inline int g_fn(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::g_rule) {
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
    return 1;
}

constexpr inline int c_fn(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::c_rule) {
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
    return 1;
}

constexpr inline int a_fn(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::a_rule) {
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
    return 1;
}

constexpr inline int u_fn(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::u_rule) {
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
    return 1;
}

constexpr inline int dg_fn(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::dg_rule) {
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
    return 1;
}

constexpr inline int dc_fn(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::dc_rule) {
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
    return 1;
}

constexpr inline int da_fn(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::da_rule) {
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
    return 1;
}

constexpr inline int dt_fn(std::string_view s2, std::string_view s3) {
    for (const auto& rule : rules::dt_rule) {
        if (rule.s2 == s2 && rule.s3 == s3) {
            return rule.result;
        }
    }
    return 1;
}

// Function pointer lookup
constexpr HandlerFunc getHandler(resTokenType type) {
  switch (type) {
    case resTokenType::PHE:
      return phe_fn;
    case resTokenType::TYR:
      return tyr_fn;
    case resTokenType::TRP:
      return trp_fn;
    // NOTE: Atom types are not unique among different force fields, so this approach is not
    // going to work as simple as this.
    case resTokenType::HIS:
    case resTokenType::HSD:
    case resTokenType::HSE:
    case resTokenType::HSP:
    case resTokenType::HID:
    case resTokenType::HIE:
    case resTokenType::HIP:
      return his_fn;
    case resTokenType::GLU:
      return glu_fn;
    case resTokenType::ASP:
      return asp_fn;
    case resTokenType::ASN:
      return asn;
    case resTokenType::GLN:
      return gln;
    case resTokenType::ARG:
      return arg;
    case resTokenType::G:
      return g_fn;
    case resTokenType::C:
      return c_fn;
    case resTokenType::A:
      return a_fn;
    case resTokenType::U:
      return u_fn;
    case resTokenType::DA:
      return da_fn;
    case resTokenType::DC:
      return dc_fn;
    case resTokenType::DG:
      return dg_fn;
    case resTokenType::DT:
      return dt_fn;
    case resTokenType::I:
    case resTokenType::N:
    case resTokenType::DI:
    case resTokenType::DU:
    case resTokenType::DN:
    case resTokenType::APN:
    case resTokenType::CPN:
    case resTokenType::TPN:
    case resTokenType::GPN:
      return default_base_fn;

    case resTokenType::SOL:
    case resTokenType::WAT:
    case resTokenType::HOH:
    case resTokenType::H2O:
    case resTokenType::W:
    case resTokenType::DOD:
    case resTokenType::D3O:
    case resTokenType::TIP:
    case resTokenType::TIP3:
    case resTokenType::TIP4:
    case resTokenType::SPC:
      return default_fn; 

    default:
      return default_amino_acid_fn;
  }
}

// Template metaprogramming to generate array at compile-time
template <std::size_t... Is>
constexpr auto generateHandlerArray(std::index_sequence<Is...>) {
  return std::array<HandlerFunc, sizeof...(Is)>{
      getHandler(static_cast<resTokenType>(Is))...};
}

// NOTE: The compile-time generation is unlikely to lead to any meaningful
// performance improvement. We should instead generate it as part of the main
// Lahuta interface as part of parsing the input file. This way the size of the
// predefined data is not a concern, and we can only generate the handlers
// (residues) that are actually needed.
// NOTE: Any performance improvements are not due to lookup times, but rather
// decreasing the number of atoms for SMARTS pattern matching.
constexpr auto handlers = generateHandlerArray(
    std::make_index_sequence<static_cast<size_t>(resTokenType::UNKNOWN)>{});

// Main processing function
inline int process(resTokenType t, std::string_view s2, std::string_view s3) {
  // `process` is responsible for swapping the arguments if necessary
  if (s2 > s3) {
    std::swap(s2, s3);
  }
  return handlers[static_cast<size_t>(t)](s2, s3);
}

struct PossiblyBonded {
  bool atom1_in_table = false;
  bool atom2_in_table = false;
  int bond_order = 0;

  PossiblyBonded(int order) : bond_order(order) {}
  PossiblyBonded(bool a1, bool a2, int order)
      : atom1_in_table(a1), atom2_in_table(a2), bond_order(order) {}

  operator int() const { return bond_order; }
  explicit operator bool() const { return bond_order != 0; }

  void setOrder(int order) { bond_order = order; }
  void setIsMetal(int atom) {
    if (atom == 1) {
      atom1_is_metal = true;
    } else {
      atom2_is_metal = true;
    }
  }

private:
  bool atom1_is_metal = false;
  bool atom2_is_metal = false;
};

inline PossiblyBonded getIntraBondOrder(RDKit::Atom *atom1, RDKit::Atom *atom2) {

  static auto is_same_conformer = [](std::string altlocA, std::string altlocB) {
    return altlocA.empty() || altlocB.empty() || altlocA == altlocB;
  };

  auto *infoA = static_cast<RDKit::AtomPDBResidueInfo *>(atom1->getMonomerInfo());
  auto *infoB = static_cast<RDKit::AtomPDBResidueInfo *>(atom2->getMonomerInfo());

  auto entryA = res_name_table(infoA->getResidueName().c_str(),
                               infoA->getResidueName().length());
  auto entryB = res_name_table(infoB->getResidueName().c_str(),
                               infoB->getResidueName().length());

  bool atom1_in_table = entryA != resTokenType::UNKNOWN;
  bool atom2_in_table = entryB != resTokenType::UNKNOWN;

  if (entryA == resTokenType::UNKNOWN || entryB == resTokenType::UNKNOWN) {
    return PossiblyBonded{atom1_in_table, atom2_in_table, 0};
  }
  if (!is_same_conformer(infoA->getAltLoc(), infoB->getAltLoc())) {
    return PossiblyBonded{atom1_in_table, atom2_in_table, 0};
  }

  int order = process(entryA, infoA->getName(), infoB->getName());
  return {atom1_in_table, atom2_in_table, order};
}

#endif // LAHUTA_BONDS_HPP
