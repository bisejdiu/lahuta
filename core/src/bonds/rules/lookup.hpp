#ifndef LAHUTA_BBONDS_HPP
#define LAHUTA_BBONDS_HPP

#include <string>
#include <string_view>

#include <rdkit/GraphMol/MonomerInfo.h>

#include "common.hpp"
#include "rules.hpp"
#include "token_lookup.hpp"

namespace lahuta {

template <typename T> struct is_rule : std::false_type {};
template <std::size_t N> struct is_rule<std::array<rules::Rule, N>> : std::true_type {};

template <auto &R> constexpr BondType rule_handler(std::string_view s2, std::string_view s3) {

  static_assert(is_rule<std::remove_cv_t<std::remove_reference_t<decltype(R)>>>::value, "Invalid rule type");

  // linear search
  for (const auto &rule : R) {
    if (rule.s2 == s2 && rule.s3 == s3) {
      return rule.result;
    }
  }

  return BondType::SINGLE;
}

using HandlerFunc = BondType (*)(std::string_view s2, std::string_view s3);
constexpr HandlerFunc get_handler(resTokenType type) {
  switch (type) {
    case resTokenType::PHE:
      return rule_handler<rules::phe_rules>;
    case resTokenType::TYR:
      return rule_handler<rules::tyr_rules>;
    case resTokenType::TRP:
      return rule_handler<rules::trp_rules>;
    // FIX: Atom types are not unique among different force fields, so this
    // approach is not going to work as simple as this.
    case resTokenType::HIS:
    case resTokenType::HSD:
    case resTokenType::HSE:
    case resTokenType::HSP:
    case resTokenType::HID:
    case resTokenType::HIE:
    case resTokenType::HIP:
      return rule_handler<rules::his_rules>;
    case resTokenType::GLU:
      return rule_handler<rules::glu_rules>;
    case resTokenType::ASP:
      return rule_handler<rules::asp_rules>;
    case resTokenType::ASN:
      return rule_handler<rules::asn_rules>;
    case resTokenType::GLN:
      return rule_handler<rules::gln_rules>;
    case resTokenType::ARG:
      return rule_handler<rules::arg_rules>;
    case resTokenType::G:
      return rule_handler<rules::g_rules>;
    case resTokenType::C:
      return rule_handler<rules::c_rules>;
    case resTokenType::A:
      return rule_handler<rules::a_rules>;
    case resTokenType::U:
      return rule_handler<rules::u_rules>;
    case resTokenType::DA:
      return rule_handler<rules::da_rules>;
    case resTokenType::DC:
      return rule_handler<rules::dc_rules>;
    case resTokenType::DG:
      return rule_handler<rules::dg_rules>;
    case resTokenType::DT:
      return rule_handler<rules::dt_rules>;
    case resTokenType::I:
    case resTokenType::N:
    case resTokenType::DI:
    case resTokenType::DU:
    case resTokenType::DN:
    case resTokenType::APN:
    case resTokenType::CPN:
    case resTokenType::TPN:
    case resTokenType::GPN:
      return rule_handler<rules::default_base_rules>;

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
      return rule_handler<rules::default_aa_rules>;
  }
}

/// generate an array of function pointers for all residue types
template <std::size_t... Is> //
constexpr auto get_handlers(std::index_sequence<Is...>) {
  return std::array<HandlerFunc, sizeof...(Is)>{get_handler(static_cast<resTokenType>(Is))...};
}

inline BondType process(resTokenType t, std::string_view s2, std::string_view s3) {
  static auto handlers = get_handlers(std::make_index_sequence<static_cast<size_t>(resTokenType::UNKNOWN)>{});
  if (s2 > s3) std::swap(s2, s3);
  return handlers[static_cast<size_t>(t)](s2, s3);
}

struct PossiblyBonded {
  bool atom1_is_predef = false;
  bool atom2_is_predef = false;
  BondType bond_type = BondType::UNSPECIFIED;
  bool is_bond_to_h = false;

  PossiblyBonded(bool a1, bool a2, BondType bt, bool is_h = false) : atom1_is_predef(a1), atom2_is_predef(a2), bond_type(bt), is_bond_to_h(is_h) {}

  explicit operator bool() const { return bond_type != BondType::UNSPECIFIED; }
};

inline PossiblyBonded getIntraBondOrder(AtomInfo a1, AtomInfo a2) {

  auto res_a = res_name_table(a1.info->getResidueName().c_str(), a1.info->getResidueName().length());
  auto res_b = res_name_table(a2.info->getResidueName().c_str(), a2.info->getResidueName().length());

  bool atom1_in_table = res_a != resTokenType::UNKNOWN;
  bool atom2_in_table = res_b != resTokenType::UNKNOWN;

  if (!atom1_in_table || !atom2_in_table) {
    return PossiblyBonded{atom1_in_table, atom2_in_table, BondType::UNSPECIFIED};
  }

  // at this point both atoms are in the table
  if (a1.is_hydrogen ^ a2.is_hydrogen) {
    return PossiblyBonded{atom1_in_table, atom2_in_table, BondType::SINGLE, true};
  }

  BondType bond_type = process(res_a, a1.info->getName(), a2.info->getName());
  return {atom1_in_table, atom2_in_table, bond_type};
}

} // namespace lahuta

#endif // LAHUTA_BBONDS_HPP
