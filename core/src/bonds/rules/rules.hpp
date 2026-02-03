/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::ostringstream os; std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   std::copy(parts.begin(), parts.end(), std::ostream_iterator<std::string_view>(os));
 *   return os.str();
 * }();
 *
 */

#ifndef LAHUTA_RULES_HPP
#define LAHUTA_RULES_HPP

#include <array>
#include <string_view>

#include <rdkit/GraphMol/Bond.h>

#include "token.h"

using BondType = RDKit::Bond::BondType;

namespace lahuta {

inline constexpr std::array<resTokenType, 28> PredefinedResidues{
    resTokenType::ALA, resTokenType::ARG, resTokenType::ASN, resTokenType::ASP, resTokenType::CYS,
    resTokenType::GLN, resTokenType::GLU, resTokenType::GLY, resTokenType::HIS, resTokenType::ILE,
    resTokenType::LEU, resTokenType::LYS, resTokenType::MET, resTokenType::PHE, resTokenType::PRO,
    resTokenType::SER, resTokenType::THR, resTokenType::TRP, resTokenType::TYR, resTokenType::VAL,
    resTokenType::A,   resTokenType::G,   resTokenType::C,   resTokenType::U,   resTokenType::DA,
    resTokenType::DC,  resTokenType::DG,  resTokenType::DT,
    /*resTokenType::HOH*/
};

namespace rules {

struct Rule {
  std::string_view s2;
  std::string_view s3;
  BondType result;
};

// clang-format off
constexpr std::array<Rule, 1> default_aa_rules = {{
    {"C", "O", BondType::DOUBLE},
}};

constexpr std::array<Rule, 1> default_base_rules = {{
    {"OP1", "P", BondType::DOUBLE},
}};

constexpr std::array<Rule, 6> his_rules = {{
    {"CE1", "NE2", BondType::AROMATIC},
    {"CD2", "NE2", BondType::AROMATIC},
    {"CE1", "ND1", BondType::AROMATIC},
    {"CD2", "CG", BondType::AROMATIC},
    {"CG", "ND1", BondType::AROMATIC},
    {"C", "O", BondType::DOUBLE}
}};

constexpr std::array<Rule, 2> arg_rules = {{
    {"CZ", "NH2", BondType::DOUBLE},
    {"C", "O", BondType::DOUBLE}
}};

constexpr std::array<Rule, 7> phe_rules = {{
    {"CD1", "CE1", BondType::AROMATIC},
    {"CD1", "CG", BondType::AROMATIC},
    {"CD2", "CE2", BondType::AROMATIC},
    {"CD2", "CG", BondType::AROMATIC},
    {"CE1", "CZ", BondType::AROMATIC},
    {"CE2", "CZ", BondType::AROMATIC},
    {"C", "O", BondType::DOUBLE}
}};

constexpr std::array<Rule, 11> trp_rules = {{
    {"CH2", "CZ3", BondType::AROMATIC},
    {"CH2", "CZ2", BondType::AROMATIC},
    {"CE3", "CZ3", BondType::AROMATIC},
    {"CD2", "CE3", BondType::AROMATIC},
    {"CD2", "CG", BondType::AROMATIC},
    {"CD2", "CE2", BondType::AROMATIC},
    {"CE2", "CZ2", BondType::AROMATIC},
    {"CE2", "NE1", BondType::AROMATIC},
    {"CD1", "NE1", BondType::AROMATIC},
    {"CD1", "CG", BondType::AROMATIC},
    {"C", "O", BondType::DOUBLE}
}};

constexpr std::array<Rule, 2> asn_rules = {{
    {"CG", "OD1", BondType::DOUBLE},
    {"C", "O", BondType::DOUBLE}
}};

constexpr std::array<Rule, 2> gln_rules = {{
    {"CD", "OE1", BondType::DOUBLE},
    {"C", "O", BondType::DOUBLE}
}};

constexpr std::array<Rule, 7> tyr_rules = {{
    {"CD1", "CE1", BondType::AROMATIC},
    {"CD1", "CG", BondType::AROMATIC},
    {"CD2", "CE2", BondType::AROMATIC},
    {"CD2", "CG", BondType::AROMATIC},
    {"CE1", "CZ", BondType::AROMATIC},
    {"CE2", "CZ", BondType::AROMATIC},
    {"C", "O", BondType::DOUBLE}
}};

constexpr std::array<Rule, 2> asp_rules = {{
    {"CG", "OD1", BondType::DOUBLE},
    {"C", "O", BondType::DOUBLE}
}};

constexpr std::array<Rule, 2> glu_rules = {{
    {"CD", "OE1", BondType::DOUBLE},
    {"C", "O", BondType::DOUBLE}
}};

constexpr std::array<Rule, 11> a_rules = {{
    {"C2", "N1", BondType::AROMATIC},
    {"C6", "N1", BondType::AROMATIC},
    {"C4", "N9", BondType::AROMATIC},
    {"C8", "N9", BondType::AROMATIC},
    {"C4", "N3", BondType::AROMATIC},
    {"C4", "C5", BondType::AROMATIC},
    {"C2", "N3", BondType::AROMATIC},
    {"C5", "C6", BondType::AROMATIC},
    {"C5", "N7", BondType::AROMATIC},
    {"C8", "N7", BondType::AROMATIC},
    {"OP1", "P", BondType::DOUBLE}
}};

constexpr std::array<Rule, 12> g_rules = {{
    {"C2", "N1", BondType::AROMATIC}, //
    {"C6", "N1", BondType::AROMATIC}, //
    {"C4", "N9", BondType::AROMATIC}, //
    {"C8", "N9", BondType::AROMATIC}, //
    {"C4", "N3", BondType::AROMATIC}, //
    {"C4", "C5", BondType::AROMATIC}, //
    {"C2", "N3", BondType::AROMATIC}, //
    {"C5", "C6", BondType::AROMATIC}, //
    {"C5", "N7", BondType::AROMATIC}, //
    {"C8", "N7", BondType::AROMATIC}, //
    {"C6", "O6", BondType::DOUBLE}, //
    {"OP1", "P", BondType::DOUBLE}
}};

constexpr std::array<Rule, 8> c_rules = {{
    {"C4", "N3", BondType::AROMATIC},
    {"C4", "C5", BondType::AROMATIC},
    {"C2", "N3", BondType::AROMATIC},
    {"C5", "C6", BondType::AROMATIC},
    {"C2", "N1", BondType::AROMATIC},
    {"C6", "N1", BondType::AROMATIC},
    {"C2", "O2", BondType::DOUBLE},
    {"OP1", "P", BondType::DOUBLE}
}};

constexpr std::array<Rule, 9> u_rules = {{
    {"C5", "C6", BondType::AROMATIC},
    {"C6", "N1", BondType::AROMATIC},
    {"C4", "C5", BondType::AROMATIC},
    {"C4", "N3", BondType::AROMATIC},
    {"C2", "N3", BondType::AROMATIC},
    {"C2", "N1", BondType::AROMATIC},
    {"C4", "O4", BondType::DOUBLE},
    {"C2", "O2", BondType::DOUBLE},
    {"OP1", "P", BondType::DOUBLE}
}};

constexpr std::array<Rule, 12> dg_rules = {{
    {"C5", "C6", BondType::AROMATIC},
    {"C6", "N1", BondType::AROMATIC},
    {"C5", "N7", BondType::AROMATIC},
    {"C4", "C5", BondType::AROMATIC},
    {"C8", "N7", BondType::AROMATIC},
    {"C8", "N9", BondType::AROMATIC},
    {"C4", "N3", BondType::AROMATIC},
    {"C4", "N9", BondType::AROMATIC},
    {"C2", "N3", BondType::AROMATIC},
    {"C2", "N1", BondType::AROMATIC},
    {"C6", "O6", BondType::DOUBLE},
    {"OP1", "P", BondType::DOUBLE}
}};

constexpr std::array<Rule, 8> dc_rules = {{
    {"C6", "N1", BondType::AROMATIC},
    {"C5", "C6", BondType::AROMATIC},
    {"C2", "N1", BondType::AROMATIC},
    {"C4", "C5", BondType::AROMATIC},
    {"C4", "N3", BondType::AROMATIC},
    {"C2", "N3", BondType::AROMATIC},
    {"C2", "O2", BondType::DOUBLE},
    {"OP1", "P", BondType::DOUBLE}
}};

constexpr std::array<Rule, 11> da_rules = {{
    {"C4", "N3", BondType::AROMATIC},
    {"C4", "C5", BondType::AROMATIC},
    {"C4", "N9", BondType::AROMATIC},
    {"C2", "N3", BondType::AROMATIC},
    {"C2", "N1", BondType::AROMATIC},
    {"C6", "N1", BondType::AROMATIC},
    {"C5", "N7", BondType::AROMATIC},
    {"C5", "C6", BondType::AROMATIC},
    {"C8", "N7", BondType::AROMATIC},
    {"C8", "N9", BondType::AROMATIC},
    {"OP1", "P", BondType::DOUBLE}
}};

constexpr std::array<Rule, 9> dt_rules = {{
    {"C4", "N3", BondType::AROMATIC},
    {"C4", "C5", BondType::AROMATIC},
    {"C2", "N3", BondType::AROMATIC},
    {"C5", "C6", BondType::AROMATIC},
    {"C2", "N1", BondType::AROMATIC},
    {"C6", "N1", BondType::AROMATIC},
    {"C2", "O2", BondType::DOUBLE},
    {"C4", "O4", BondType::DOUBLE},
    {"OP1", "P", BondType::DOUBLE}
}};

} // namespace rules
} // namespace lahuta
//
#endif // LAHUTA_RULES_HPP
