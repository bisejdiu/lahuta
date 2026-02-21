/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto f = [](auto arg) {
 *     using T = std::decay_t<decltype(arg)>;
 *     if constexpr (std::disjunction_v<std::is_same<T, const char*>, std::is_same<T, std::string_view>>) return std::string(arg);
 *   };
 *   return f("besian") + f("sejdiu") + f("@gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_TYPING_SMARTS_MATCHING_HPP
#define LAHUTA_TYPING_SMARTS_MATCHING_HPP

#include "typing/types.hpp"

namespace lahuta {

constexpr std::pair<const char *, AtomType> AtomTypeSMARTS[] = {
    {"[$([nH]:@c(=O))]", 1029_at},
    {"[$([N;H2;v3;$(N-C(=O))])]", 1029_at},
    {"[$([n;H1;v3;!$([nH]cccc)])]", 1029_at},
    // Requires RingInfo initialization
    {"[#8,#9,$([#16;H0,H1;v2,v1]),$([N;v3;!$(N-*"
     "=!@[O,N,P,S]);!$(N-!@a);!$([NH]=!@*)]),$([nH0;+0])]",
     1029_at},
    /**/
    {"[$(n:a:[nH])]", 2_at},
    {"[$([O;H0;$(O=C-[NH2])])]", 2_at},
    {"[$([O;H0;$(O=C([OH])-*)])]", 2_at},
    {"[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]", 2_at},

    {"[Cl,Br,I;X1;$([Cl,Br,I]-[#6])]", 10_at},

    {"[$([N;H2&+0][C;!$(C=*)]),$([N;H1&+0]([C;!$(C=*)])[C;!$(C=*)]),$([N;H0&+0]"
     "([C;!$(C=*)])([C;!$(C=*)])[C;!$(C=*)]);!$(N[a])]",
     16_at},
    {"NC(=N)", 16_at},
    {"[#7;+;!$([N+]-[O-])]", 16_at},
    {"[$([*+1,*+2,*+3]);!$([N+]-[O-])]", 16_at},

    {"[*-1,*-2]", 32_at},
    {"[$([OH,O-]-[C,S,N,P,Cl,Br,I]=O),$(O=[C,S,N,P,Cl,Br,I]-[OH,O-])]", 32_at},

    {"[$([OH0]=[CX3,c]);!$([OH0]=[CX3,c]-[OH,O-])]", 64_at},
    {"[$([CX3,c]=[OH0]);!$([CX3,c](=[OH0])-[OH,O-])]", 128_at},
    {"[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,Cl+0,Br+0,I+0]", 512_at},

    // Require RingInfo initialization
    {"[n;R1]1[c;R1][n;R1][c;R1][c;R1]1", 16_at},
    {"[Xx]", 1_at},

    {"[#6!H0]", 8_at},
    {"[Cl,Br,I;X1;$([Cl,Br,I]-[#6])]", 8_at},
};

inline std::vector<AtomType> match_atom_types(RDKit::ROMol &mol) {
  thread_local std::array<RDKit::ROMol*, std::size(AtomTypeSMARTS)> patterns = [] {
    std::array<RDKit::ROMol*, std::size(AtomTypeSMARTS)> temp{};
    for (size_t i = 0; i < std::size(AtomTypeSMARTS); ++i) {
      temp[i] = RDKit::SmartsToMol(AtomTypeSMARTS[i].first);
    }
    return temp;
  }();

  RDKit::SubstructMatchParameters params;
  params.maxMatches = mol.getNumAtoms();

  std::vector<AtomType> types(mol.getNumAtoms(), AtomType::None);
  for (size_t i = 0; i < std::size(AtomTypeSMARTS); ++i) {
    const auto &[smarts, atom_type] = AtomTypeSMARTS[i];
    RDKit::ROMol *pattern = patterns[i];

    SubStrMatches match_list;
    RDKit::SubstructMatch(mol, *pattern, match_list);

    for (const auto &match : match_list) {
      for (const auto &pair : match) {
        types[pair.second] |= atom_type;
      }
    }
  }

  return types;
}

} // namespace lahuta

#endif // LAHUTA_TYPING_SMARTS_MATCHING_HPP
