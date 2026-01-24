#ifndef LAHUTA_DEFINITIONS_HPP
#define LAHUTA_DEFINITIONS_HPP

#include <array>
#include <string>
#include <unordered_set>
#include <vector>

#include "bonds/rules/rules.hpp"
#include "bonds/rules/token_lookup.hpp"

namespace std {
template <> struct hash<lahuta::resTokenType> {
  size_t operator()(const lahuta::resTokenType e) const noexcept { return static_cast<size_t>(e); }
};
} // namespace std

namespace lahuta {
namespace definitions {

template <size_t N>
inline std::function<bool(const std::string &)> //
make_tester(const std::array<resTokenType, N> &values) {
  return [values](const std::string &res_name) {
    auto entry = res_name_table(res_name.c_str(), res_name.length());
    return std::find(values.begin(), values.end(), entry) != values.end();
  };
}

inline std::function<bool(const std::string &)>
make_tester(resTokenType range_start, resTokenType range_end) {
  return [range_start, range_end](const std::string &res_name) {
    auto entry = res_name_table(res_name.c_str(), res_name.length());
    return entry >= range_start && entry <= range_end;
  };
}

template <size_t N>
inline std::function<bool(const std::string &)>
make_tester(const std::array<std::pair<resTokenType, resTokenType>, N> &ranges) {
  return [ranges](const std::string &res_name) {
    auto entry = res_name_table(res_name.c_str(), res_name.length());
    for (const auto &[start, end] : ranges) {
      if (entry >= start && entry <= end) {
        return true;
      }
    }
    return false;
  };
}

// clang-format off
inline constexpr std::array<resTokenType, 3> PositiveChargedResidues = {
    resTokenType::ARG, resTokenType::HIS, resTokenType::LYS
};

inline constexpr std::array<resTokenType, 2> NegativeChargedResidues = {
    resTokenType::GLU, resTokenType::ASP
};

inline constexpr std::array<std::pair<resTokenType, resTokenType>, 2> HistidineResidues = {
    std::make_pair(resTokenType::HIS, resTokenType::HIS),
    std::make_pair(resTokenType::HSD, resTokenType::HIP)
};

inline constexpr std::pair<resTokenType, resTokenType> StandardAminoAcids = {
    resTokenType::GLY, resTokenType::HYP
};

inline constexpr std::pair<resTokenType, resTokenType> PolymerResiduesRange = {
    resTokenType::GLY, resTokenType::GPN
};

inline constexpr std::pair<resTokenType, resTokenType> BaseResiduesRange = {
    resTokenType::GLY, resTokenType::GPN
};

inline constexpr int _AROMATIC_RESIDUES_COUNT = 28;
inline constexpr std::array<resTokenType, _AROMATIC_RESIDUES_COUNT> AromaticResidues = {
    resTokenType::PHE, resTokenType::TYR, // phenylalanine, tyrosine
    resTokenType::HIS, resTokenType::TRP, // histidine, tryptophan
    resTokenType::DHI, resTokenType::DPN, // d-histidine, d-phenylalanine
    resTokenType::DTR, resTokenType::DTY, // d-tryptophan, d-tyrosine
    resTokenType::HSD, resTokenType::HSE, // histidine delta, epsilon
    resTokenType::HSP, resTokenType::HID, // histidine protonated, delta
    resTokenType::HIE, resTokenType::HIP, // histidine epsilon, protonated
    resTokenType::PTR, resTokenType::A,   // phosphotyrosine, adenine
    resTokenType::G,   resTokenType::C,   // guanine, cytosine
    resTokenType::U,   resTokenType::I,   // uracil, inosine
    resTokenType::N,   resTokenType::DA,  // nucleotide, deoxyadenosine
    resTokenType::DC,  resTokenType::DG,  // deoxycytidine, deoxyguanosine
    resTokenType::DT,  resTokenType::DI,  // deoxythymidine, deoxyinosine
    resTokenType::DU,  resTokenType::DN   // deoxyuridine, deoxynucleotide
};


namespace arom_rings {
// old-style enum to allow implicit conversion to int.
enum RingSize {
    RS_None = 0,
    RS_3    = 1 << 0,
    RS_4    = 1 << 1,
    RS_5    = 1 << 2,
    RS_6    = 1 << 3,
    RS_7    = 1 << 4,
    RS_8    = 1 << 5,
};

//
// Using an array forces a linear search, but: 1. allows for optimizations (short-circuiting) and
// (2) forces the same size with _AromaticResidues_ (which is important) and (3) for the current size it
// actually should be comparable (likely even faster) than a hash map.  -Besian, March 2025
//
constexpr std::array<std::pair<const char *, RingSize>, _AROMATIC_RESIDUES_COUNT> AromaticResiduesRingSizes {{
  {"PHE", RS_6},
  {"TYR", RS_6},
  {"HIS", RS_5},
  {"TRP", static_cast<RingSize>(RS_5 | RS_6)},

  {"A",   static_cast<RingSize>(RS_5 | RS_6)},
  {"G",   RS_6},
  {"C",   RS_6},
  {"U",   RS_6},

  {"DA",  RS_6},
  {"DC",  RS_6},
  {"DG",  RS_6},
  {"DT",  RS_6},

  {"N",   RS_6},
  {"I",   RS_6},

  {"DN",  RS_6},
  {"DU",  RS_6},
  {"DI",  RS_6},

  {"DTR", static_cast<RingSize>(RS_5 | RS_6)},
  {"DTY", static_cast<RingSize>(RS_5 | RS_6)},
  {"DHI", RS_5},
  {"DPN", RS_5},

  {"HSD", RS_5},
  {"HSE", RS_5},
  {"HSP", RS_5},
  {"HID", RS_5},
  {"HIE", RS_5},
  {"HIP", RS_5},
  {"PTR", RS_6}
}};

inline std::vector<int> get_ringsizes(RingSize ring_mask) {
      std::vector<int> sizes;
    if (ring_mask & RS_3) sizes.push_back(3);
    if (ring_mask & RS_4) sizes.push_back(4);
    if (ring_mask & RS_5) sizes.push_back(5);
    if (ring_mask & RS_6) sizes.push_back(6);
    if (ring_mask & RS_7) sizes.push_back(7);
    if (ring_mask & RS_8) sizes.push_back(8);
    return sizes;
}

} // namespace arom_rings


const std::unordered_set<std::string> ProteinBackboneAtoms = {
  "CA",  "C",   "N",   "O",   "O1",  "O2", "OC1",
  "OC2", "OT1", "OT2", "OX1", "OXT", "H",  "H1",
  "H2",  "H3",  "HA",  "HN",  "HXT", "BB"
};

const std::unordered_set<std::string> NucleicBackboneAtoms = {
  "P", "OP1", "OP2", "HOP2", "HOP3", "O2\'", "O3\'",
  "O4\'", "O5\'", "C1\'", "C2\'", "C3\'", "C4\'",
  "C5\'", "H1\'", "H2\'", "H2\'\'", "HO2\'", "H3\'",
  "H4\'", "H5\'", "H5\'\'", "HO3\'", "HO5\'", "O2*",
  "O3*", "O4*", "O5*", "C1*", "C2*", "C3*", "C4*", "C5*"
};

/*using ResTesterFunc = std::function<bool(const std::string &)>;*/
const auto is_water   = make_tester(resTokenType::SOL, resTokenType::SPC);
const auto is_rna     = make_tester(resTokenType::A, resTokenType::N);
const auto is_dna     = make_tester(resTokenType::DA, resTokenType::DN);
const auto is_base    = make_tester(BaseResiduesRange.first, BaseResiduesRange.second);
const auto is_nucleic = make_tester(resTokenType::A, resTokenType::DN);
const auto is_polymer = make_tester(PolymerResiduesRange.first, PolymerResiduesRange.second);

const auto is_histidine        = make_tester(HistidineResidues);
const auto is_positive_charge  = make_tester(PositiveChargedResidues);
const auto is_negative_charge  = make_tester(NegativeChargedResidues);
const auto is_standard_protein = make_tester(resTokenType::GLY, resTokenType::ARG);
const auto is_protein_extended = make_tester(StandardAminoAcids.first, StandardAminoAcids.second);

const auto is_aromatic   = make_tester(AromaticResidues);
const auto is_predefined = make_tester(PredefinedResidues);

} // namespace definitions
} // namespace lahuta

#endif // LAHUTA_DEFINITIONS_HPP
