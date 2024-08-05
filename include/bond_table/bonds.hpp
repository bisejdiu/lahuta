#ifndef LAHUTA_BONDS_HPP
#define LAHUTA_BONDS_HPP
#include <string>

#include "token-gperf-generated.hpp"
#include "token.h"

#include "GraphMol/Atom.h"
#include "GraphMol/MonomerInfo.h"

using namespace lahuta;

// Function pointer type for our handlers
using HandlerFunc = int (*)(std::string_view s2, std::string_view s3);

inline int defaultHandler(std::string_view s2, std::string_view s3) {
  return 1;
}
inline int defaultAminoAcidHandler(std::string_view s2, std::string_view s3) {
  if (s3 == "O" && s2 == "C") {
    return 2;
  }
  return 1;
}
inline int defaultBaseHandler(std::string_view s2, std::string_view s3) {
  if (s2 == "OP1" && s3 == "P") {
    return 2;
  }
  return 1;
}

// Implementation of handler functions
inline int HIS(std::string_view s2, std::string_view s3) {
  if (s2 == "CD2" && s3 == "CG") {
    return 2;
  }
  if (s2 == "CE1" && s3 == "ND1") {
    return 2;
  }
  if (s3 == "O" && s2 == "C") {
    return 2;
  }
  return 1;
}

inline int ARG(std::string_view s2, std::string_view s3) {
  if (s2 == "CZ" && s3 == "NH2") {
    return 2;
  }
  if (s3 == "O" && s2 == "C") {
    return 2;
  }
  return 1;
}

inline int PHE(std::string_view s2, std::string_view s3) {
  if (s2 == "CE1" && s3 == "CZ") {
    return 2;
  }
  if (s2 == "CD2" && s3 == "CE2") {
    return 2;
  }
  if (s2 == "CD1" && s3 == "CG") {
    return 2;
  }
  if (s3 == "O" && s2 == "C") {
    return 2;
  }
  return 1;
}

inline int TRP(std::string_view s2, std::string_view s3) {
  if (s2 == "CD1" && s3 == "CG") {
    return 2;
  }
  if (s2 == "CD2" && s3 == "CE2") {
    return 2;
  }
  if (s2 == "CE3" && s3 == "CZ3") {
    return 2;
  }
  if (s2 == "CH2" && s3 == "CZ2") {
    return 2;
  }
  if (s3 == "O" && s2 == "C") {
    return 2;
  }
  return 1;
}

inline int ASN(std::string_view s2, std::string_view s3) {
  if (s2 == "CG" && s3 == "OD1") {
    return 2;
  }
  if (s3 == "O" && s2 == "C") {
    return 2;
  }
  return 1;
}

inline int GLN(std::string_view s2, std::string_view s3) {
  if (s2 == "CD" && s3 == "OE1") {
    return 2;
  }
  if (s3 == "O" && s2 == "C") {
    return 2;
  }
  return 1;
}

inline int TYR(std::string_view s2, std::string_view s3) {
  if (s2 == "CD1" && s3 == "CG") {
    return 2;
  }
  if (s2 == "CD2" && s3 == "CE2") {
    return 2;
  }
  if (s2 == "CE1" && s3 == "CZ") {
    return 2;
  }
  if (s3 == "O" && s2 == "C") {
    return 2;
  }
  return 1;
}

inline int ASP(std::string_view s2, std::string_view s3) {
  if (s2 == "CG" && s3 == "OD1") {
    return 2;
  }
  if (s3 == "O" && s2 == "C") {
    return 2;
  }
  return 1;
}

inline int GLU(std::string_view s2, std::string_view s3) {
  if (s2 == "CD" && s3 == "OE1") {
    return 2;
  }
  if (s3 == "O" && s2 == "C") {
    return 2;
  }
  return 1;
}

inline int G(std::string_view s2, std::string_view s3) {
  if (s2 == "C8" && s3 == "N7") {
    return 2;
  }
  if (s2 == "C4" && s3 == "C5") {
    return 2;
  }
  if (s2 == "C2" && s3 == "N3") {
    return 2;
  }
  if (s2 == "C6" && s3 == "O6") {
    return 2;
  }
  if (s2 == "OP1" && s3 == "P") {
    return 2;
  }
  return 1;
}

inline int C(std::string_view s2, std::string_view s3) {
  if (s2 == "C4" && s3 == "N3") {
    return 2;
  }
  if (s2 == "C5" && s3 == "C6") {
    return 2;
  }
  if (s2 == "C2" && s3 == "O2") {
    return 2;
  }
  if (s2 == "OP1" && s3 == "P") {
    return 2;
  }
  return 1;
}

inline int A(std::string_view s2, std::string_view s3) {
  if (s2 == "C2" && s3 == "N3") {
    return 2;
  }
  if (s2 == "C6" && s3 == "N1") {
    return 2;
  }
  if (s2 == "C4" && s3 == "C5") {
    return 2;
  }
  if (s2 == "C8" && s3 == "N7") {
    return 2;
  }
  if (s2 == "OP1" && s3 == "P") {
    return 2;
  }
  return 1;
}

inline int U(std::string_view s2, std::string_view s3) {
  if (s2 == "C5" && s3 == "C6") {
    return 2;
  }
  if (s2 == "C2" && s3 == "O2") {
    return 2;
  }
  if (s2 == "C4" && s3 == "O4") {
    return 2;
  }
  if (s2 == "OP1" && s3 == "P") {
    return 2;
  }
  return 1;
}

inline int DG(std::string_view s2, std::string_view s3) {
  if (s2 == "C8" && s3 == "N7") {
    return 2;
  }
  if (s2 == "C4" && s3 == "C5") {
    return 2;
  }
  if (s2 == "C2" && s3 == "N3") {
    return 2;
  }
  if (s2 == "C6" && s3 == "O6") {
    return 2;
  }
  if (s2 == "OP1" && s3 == "P") {
    return 2;
  }
  return 1;
}

inline int DC(std::string_view s2, std::string_view s3) {
  if (s2 == "C4" && s3 == "N3") {
    return 2;
  }
  if (s2 == "C5" && s3 == "C6") {
    return 2;
  }
  if (s2 == "C2" && s3 == "O2") {
    return 2;
  }
  if (s2 == "OP1" && s3 == "P") {
    return 2;
  }
  return 1;
}

inline int DA(std::string_view s2, std::string_view s3) {
  if (s2 == "C2" && s3 == "N3") {
    return 2;
  }
  if (s2 == "C6" && s3 == "N1") {
    return 2;
  }
  if (s2 == "C4" && s3 == "C5") {
    return 2;
  }
  if (s2 == "C8" && s3 == "N7") {
    return 2;
  }
  if (s2 == "OP1" && s3 == "P") {
    return 2;
  }
  return 1;
}

inline int DT(std::string_view s2, std::string_view s3) {
  if (s2 == "C5" && s3 == "C6") {
    return 2;
  }
  if (s2 == "C2" && s3 == "O2") {
    return 2;
  }
  if (s2 == "C4" && s3 == "O4") {
    return 2;
  }
  if (s2 == "OP1" && s3 == "P") {
    return 2;
  }
  return 1;
}

constexpr HandlerFunc getHandler(resTokenType type) {
  switch (type) {
  case resTokenType::PHE:
    return PHE;
  case resTokenType::TYR:
    return TYR;
  case resTokenType::TRP:
    return TRP;
  // NOTE: Atom types are not unique among different force fields, so this approach is not
  // going to work as simple as this.
  case resTokenType::HIS:
  case resTokenType::HSD:
  case resTokenType::HSE:
  case resTokenType::HSP:
  case resTokenType::HID:
  case resTokenType::HIE:
  case resTokenType::HIP:
    return HIS;
  case resTokenType::GLU:
    return GLU;
  case resTokenType::ASP:
    return ASP;
  case resTokenType::ASN:
    return ASN;
  case resTokenType::GLN:
    return GLN;
  case resTokenType::ARG:
    return ARG;
  case resTokenType::G:
    return G;
  case resTokenType::C:
    return C;
  case resTokenType::A:
    return A;
  case resTokenType::U:
    return U;
  case resTokenType::DA:
    return DA;
  case resTokenType::DC:
    return DC;
  case resTokenType::DG:
    return DG;
  case resTokenType::DT:
    return DT;
  case resTokenType::I:
  case resTokenType::N:
  case resTokenType::DI:
  case resTokenType::DU:
  case resTokenType::DN:
  case resTokenType::APN:
  case resTokenType::CPN:
  case resTokenType::TPN:
  case resTokenType::GPN:
    return defaultBaseHandler;

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
    return defaultHandler;

  default:
    return defaultAminoAcidHandler;
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
  // now `process` is responsible for swapping the arguments if necessary
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

inline bool is_same_conformer(std::string altlocA, std::string altlocB) {
  return altlocA.empty() || altlocB.empty() || altlocA == altlocB;
}

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
    return PossiblyBonded{atom1_in_table, atom2_in_table, 0};
  }
  if (!is_same_conformer(infoA->getAltLoc(), infoB->getAltLoc())) {
    return PossiblyBonded{atom1_in_table, atom2_in_table, 0};
  }

  int order = process(entryA, infoA->getName(), infoB->getName());
  return {atom1_in_table, atom2_in_table, order};
}

#endif // LAHUTA_BONDS_HPP
