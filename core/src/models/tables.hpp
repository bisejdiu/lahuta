/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::ostringstream os; os << "besian" << "sejdiu" << "@gmail.com";
 *   return os.str();
 * }();
 *
 */

#ifndef LAHUTA_MODEL_TABLES_HPP
#define LAHUTA_MODEL_TABLES_HPP

#include "constants.hpp"

// clang-format off
namespace lahuta {

struct AminoAcidEntry {
  const char *const *atoms;
  size_t size;
  const char *name;
  const int *ih;
  const int *at;
  const int *hyb;

  template <size_t N>
  constexpr AminoAcidEntry(
      const std::array<const char *, N> &arr,
      const char *n,
      const std::array<int, N> &ih_arr,
      const std::array<int, N> &at_arr,
      const std::array<int, N> &hyb_arr = {})
      : atoms(arr.data()), size(N), name(n), ih(ih_arr.data()), at(at_arr.data()), hyb(hyb_arr.data()) {}

  constexpr AminoAcidEntry() : atoms(nullptr), size(0), name(nullptr), ih(nullptr), at(nullptr), hyb(nullptr) {}
};


class AminoAcidTable {
private:
  static constexpr char BASE_CHAR = 'A';
  static constexpr size_t MAX_AMINO_ACID = 'Y';                        // ASCII 89
  static constexpr size_t TABLE_SIZE = MAX_AMINO_ACID - BASE_CHAR + 1; // 25

  std::array<AminoAcidEntry, TABLE_SIZE> data{};

  // just in case, we don't actually check for validity
    static constexpr std::array<char, 20> valid_codes = {
        'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
        'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
    };

public:
  constexpr AminoAcidTable() {
    data['G' - BASE_CHAR] = AminoAcidEntry(gly, "GLY", gly_ih, gly_at, gly_hyb);
    data['A' - BASE_CHAR] = AminoAcidEntry(ala, "ALA", ala_ih, ala_at, ala_hyb);
    data['V' - BASE_CHAR] = AminoAcidEntry(val, "VAL", val_ih, val_at, val_hyb);
    data['L' - BASE_CHAR] = AminoAcidEntry(leu, "LEU", leu_ih, leu_at, leu_hyb);
    data['I' - BASE_CHAR] = AminoAcidEntry(ile, "ILE", ile_ih, ile_at, ile_hyb);
    data['S' - BASE_CHAR] = AminoAcidEntry(ser, "SER", ser_ih, ser_at, ser_hyb);
    data['T' - BASE_CHAR] = AminoAcidEntry(thr, "THR", thr_ih, thr_at, thr_hyb);
    data['C' - BASE_CHAR] = AminoAcidEntry(cys, "CYS", cys_ih, cys_at, cys_hyb);
    data['M' - BASE_CHAR] = AminoAcidEntry(met, "MET", met_ih, met_at, met_hyb);
    data['P' - BASE_CHAR] = AminoAcidEntry(pro, "PRO", pro_ih, pro_at, pro_hyb);
    data['F' - BASE_CHAR] = AminoAcidEntry(phe, "PHE", phe_ih, phe_at, phe_hyb);
    data['Y' - BASE_CHAR] = AminoAcidEntry(tyr, "TYR", tyr_ih, tyr_at, tyr_hyb);
    data['W' - BASE_CHAR] = AminoAcidEntry(trp, "TRP", trp_ih, trp_at, trp_hyb);
    data['H' - BASE_CHAR] = AminoAcidEntry(his, "HIS", his_ih, his_at, his_hyb);
    data['E' - BASE_CHAR] = AminoAcidEntry(glu, "GLU", glu_ih, glu_at, glu_hyb);
    data['D' - BASE_CHAR] = AminoAcidEntry(asp, "ASP", asp_ih, asp_at, asp_hyb);
    data['N' - BASE_CHAR] = AminoAcidEntry(asn, "ASN", asn_ih, asn_at, asn_hyb);
    data['Q' - BASE_CHAR] = AminoAcidEntry(gln, "GLN", gln_ih, gln_at, gln_hyb);
    data['K' - BASE_CHAR] = AminoAcidEntry(lys, "LYS", lys_ih, lys_at, lys_hyb);
    data['R' - BASE_CHAR] = AminoAcidEntry(arg, "ARG", arg_ih, arg_at, arg_hyb);
  }

  constexpr const AminoAcidEntry &operator[](char c) const {
    // no bounds checking
    // if (c < BASE_CHAR || c > MAX_AMINO_ACID) {
    //     throw std::out_of_range("Invalid amino acid code");
    // }

    const auto &entry = data[c - BASE_CHAR];
    // if (entry.size == 0) throw std::out_of_range("Invalid amino acid code");

    return entry;
  }

  // if we need bounds checking we can use this
  constexpr bool is_valid(char c) const {
    if (c < BASE_CHAR || c > MAX_AMINO_ACID) {
      return false;
    }
    return data[c - BASE_CHAR].size > 0;
  }

  static constexpr bool is_aromatic(const char *c) {
    return c[0] == 'F' || c[0] == 'Y' || c[0] == 'W' || c[0] == 'H';
  }

  static constexpr bool is_trp(const char *c) { return c[0] == 'W'; }
};


struct AminoAcidEdges {
  const Edge *edges;
  size_t size;

  template <size_t N> constexpr AminoAcidEdges(const std::array<Edge, N> &arr) : edges(arr.data()), size(N) {}

  constexpr AminoAcidEdges() : edges(nullptr), size(0) {}
};

class AminoAcidEdgeTable {
private:
  static constexpr char BASE_CHAR = 'A';
  static constexpr size_t MAX_AMINO_ACID = 'Y';
  static constexpr size_t TABLE_SIZE = MAX_AMINO_ACID - BASE_CHAR + 1;

  std::array<AminoAcidEdges, TABLE_SIZE> data{};

public:
  constexpr AminoAcidEdgeTable() {
    data['G' - BASE_CHAR] = AminoAcidEdges(gly_edges);
    data['A' - BASE_CHAR] = AminoAcidEdges(ala_edges);
    data['V' - BASE_CHAR] = AminoAcidEdges(val_edges);
    data['L' - BASE_CHAR] = AminoAcidEdges(leu_edges);
    data['I' - BASE_CHAR] = AminoAcidEdges(ile_edges);
    data['S' - BASE_CHAR] = AminoAcidEdges(ser_edges);
    data['T' - BASE_CHAR] = AminoAcidEdges(thr_edges);
    data['C' - BASE_CHAR] = AminoAcidEdges(cys_edges);
    data['M' - BASE_CHAR] = AminoAcidEdges(met_edges);
    data['P' - BASE_CHAR] = AminoAcidEdges(pro_edges);
    data['F' - BASE_CHAR] = AminoAcidEdges(phe_edges);
    data['Y' - BASE_CHAR] = AminoAcidEdges(tyr_edges);
    data['W' - BASE_CHAR] = AminoAcidEdges(trp_edges);
    data['H' - BASE_CHAR] = AminoAcidEdges(his_edges);
    data['E' - BASE_CHAR] = AminoAcidEdges(glu_edges);
    data['D' - BASE_CHAR] = AminoAcidEdges(asp_edges);
    data['N' - BASE_CHAR] = AminoAcidEdges(asn_edges);
    data['Q' - BASE_CHAR] = AminoAcidEdges(gln_edges);
    data['K' - BASE_CHAR] = AminoAcidEdges(lys_edges);
    data['R' - BASE_CHAR] = AminoAcidEdges(arg_edges);
  }

  constexpr const AminoAcidEdges &operator[](char c) const {
    // no bounds checking
    // if (c < BASE_CHAR || c > MAX_AMINO_ACID) {
    //     throw std::out_of_range("Invalid amino acid code");
    // }
    return data[c - BASE_CHAR];
  }
};

constexpr AminoAcidEdgeTable StandardAminoAcidBondTable;
constexpr AminoAcidTable     StandardAminoAcidDataTable;

} // namespace lahuta

#endif // LAHUTA_MODEL_TABLES_HPP
