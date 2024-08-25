#ifndef LAHUTA_BBONDS_HPP
#define LAHUTA_BBONDS_HPP

#include <iterator>
#include <string>
#include <string_view>

#include "GraphMol/Bond.h"
#include "bonds/token.h"
#include "rules.hpp"
#include "token-gperf-generated.hpp"

#include "GraphMol/Atom.h"
#include "GraphMol/MonomerInfo.h"

namespace lahuta {

// using namespace RDKit;
typedef RDKit::Bond::BondType BondType;

// Function pointer type for bond rules
using HandlerFunc = BondType (*)(std::string_view s2, std::string_view s3);

template <typename T> struct is_rule_array : std::false_type {};

template <std::size_t N>
struct is_rule_array<std::array<rules::Rule, N>> : std::true_type {};

template <const auto &T>
constexpr BondType ruleHandler(std::string_view s2, std::string_view s3) {

  static_assert(
      is_rule_array<
          std::remove_cv_t<std::remove_reference_t<decltype(T)>>>::value,
      "RuleSet must be a std::array of rules::Rule");

  for (const auto &rule : T) {
    if (rule.s2 == s2 && rule.s3 == s3) {
      return rule.result;
    }
  }
  return BondType::SINGLE;
}

// Function pointer lookup
constexpr HandlerFunc getHandler(resTokenType type) {
  switch (type) {
  case resTokenType::PHE:
    return ruleHandler<rules::phe_rules>;
  case resTokenType::TYR:
    return ruleHandler<rules::tyr_rules>;
  case resTokenType::TRP:
    return ruleHandler<rules::trp_rules>;
  // NOTE: Atom types are not unique among different force fields, so this
  // approach is not going to work as simple as this.
  case resTokenType::HIS:
  case resTokenType::HSD:
  case resTokenType::HSE:
  case resTokenType::HSP:
  case resTokenType::HID:
  case resTokenType::HIE:
  case resTokenType::HIP:
    return ruleHandler<rules::his_rules>;
  case resTokenType::GLU:
    return ruleHandler<rules::glu_rules>;
  case resTokenType::ASP:
    return ruleHandler<rules::asp_rules>;
  case resTokenType::ASN:
    return ruleHandler<rules::asn_rules>;
  case resTokenType::GLN:
    return ruleHandler<rules::gln_rules>;
  case resTokenType::ARG:
    return ruleHandler<rules::arg_rules>;
  case resTokenType::G:
    return ruleHandler<rules::g_rules>;
  case resTokenType::C:
    return ruleHandler<rules::c_rules>;
  case resTokenType::A:
    return ruleHandler<rules::a_rules>;
  case resTokenType::U:
    return ruleHandler<rules::u_rules>;
  case resTokenType::DA:
    return ruleHandler<rules::da_rules>;
  case resTokenType::DC:
    return ruleHandler<rules::dc_rules>;
  case resTokenType::DG:
    return ruleHandler<rules::dg_rules>;
  case resTokenType::DT:
    return ruleHandler<rules::dt_rules>;
  case resTokenType::I:
  case resTokenType::N:
  case resTokenType::DI:
  case resTokenType::DU:
  case resTokenType::DN:
  case resTokenType::APN:
  case resTokenType::CPN:
  case resTokenType::TPN:
  case resTokenType::GPN:
    return ruleHandler<rules::default_base_rules>;

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
    return [](std::string_view, std::string_view) { return BondType::SINGLE; };

  default:
    return ruleHandler<rules::default_aa_rules>;
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
inline BondType process(resTokenType t, std::string_view s2,
                        std::string_view s3) {
  // `process` is responsible for swapping the arguments if necessary
  if (s2 > s3) {
    std::swap(s2, s3);
  }
  return handlers[static_cast<size_t>(t)](s2, s3);
}

struct PossiblyBonded {
  bool atom1_is_predef = false;
  bool atom2_is_predef = false;
  BondType bond_type = BondType::UNSPECIFIED;

  PossiblyBonded() = default;
  PossiblyBonded(BondType bt) : bond_type(bt){};
  PossiblyBonded(bool a1, bool a2, BondType bt)
      : atom1_is_predef(a1), atom2_is_predef(a2), bond_type(bt) {}

  explicit operator bool() const { return bond_type != BondType::UNSPECIFIED; }
};

inline auto is_same_conformer(const std::string &altlocA,
                              const std::string &altlocB) {
  return altlocA.empty() || altlocB.empty() || altlocA == altlocB;
};

inline PossiblyBonded getIntraBondOrder(RDKit::Atom *atom1,
                                        RDKit::Atom *atom2) {

  auto *infoA =
      static_cast<RDKit::AtomPDBResidueInfo *>(atom1->getMonomerInfo());
  auto *infoB =
      static_cast<RDKit::AtomPDBResidueInfo *>(atom2->getMonomerInfo());

  auto entryA = res_name_table(infoA->getResidueName().c_str(),
                               infoA->getResidueName().length());
  auto entryB = res_name_table(infoB->getResidueName().c_str(),
                               infoB->getResidueName().length());

  bool atom1_in_table = entryA != resTokenType::UNKNOWN;
  bool atom2_in_table = entryB != resTokenType::UNKNOWN;

  if (entryA == resTokenType::UNKNOWN || entryB == resTokenType::UNKNOWN) {
    return PossiblyBonded{atom1_in_table, atom2_in_table,
                          BondType::UNSPECIFIED};
  }
  if (!is_same_conformer(infoA->getAltLoc(), infoB->getAltLoc())) {
    return PossiblyBonded{atom1_in_table, atom2_in_table,
                          BondType::UNSPECIFIED};
  }

  BondType bond_type = process(entryA, infoA->getName(), infoB->getName());
  PossiblyBonded pb{atom1_in_table, atom2_in_table, bond_type};

  return pb;
}

} // namespace lahuta

#endif // LAHUTA_BBONDS_HPP
