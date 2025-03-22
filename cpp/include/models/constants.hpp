#ifndef LAHUTA_MODEL_CONSTANTS_HPP
#define LAHUTA_MODEL_CONSTANTS_HPP

#include "GraphMol/Bond.h"
#include <array>

using BondType = RDKit::Bond::BondType;

namespace lahuta {

//
//  The first N has 3 implicit Hs. Others have 1.
//  CA has 1 implicit H, except glycine which has 2.
//

constexpr std::array<const char*, 4>  gly = {"N", "CA", "C", "O"};
constexpr std::array<int, 4>  gly_ih = {1, 2, 0, 0};
constexpr std::array<int, 4>  gly_at = {1026, 8, 0, 9217};

constexpr std::array<const char*, 5>  ala = {"N", "CA", "C", "CB", "O"};
constexpr std::array<int, 5>  ala_ih = {1, 1, 0, 3, 0};
constexpr std::array<int, 5>  ala_at = {1026, 8, 0, 512, 9217};

constexpr std::array<const char*, 7>  val = {"N", "CA", "C", "CB", "O", "CG1", "CG2"};
constexpr std::array<int, 7>  val_ih = {1, 1, 0, 1, 0, 3, 3};
constexpr std::array<int, 7>  val_at = {1026, 8, 0, 512, 9217, 512, 512};

constexpr std::array<const char*, 8>  leu = {"N", "CA", "C", "CB", "O", "CG", "CD1", "CD2"};
constexpr std::array<int, 8>  leu_ih = {1, 1, 0, 2, 0, 1, 3, 3};
constexpr std::array<int, 8>  leu_at = {1026, 8, 0, 512, 9217, 512, 512, 512};

constexpr std::array<const char*, 8>  ile = {"N", "CA", "C", "CB", "O", "CG1", "CG2", "CD1"};
constexpr std::array<int, 8>  ile_ih = {1, 1, 0, 1, 0, 2, 3, 3};
constexpr std::array<int, 8>  ile_at = {1026, 8, 0, 512, 9217, 512, 512, 512};

constexpr std::array<const char*, 6>  ser = {"N", "CA", "C", "CB", "O", "OG"};
constexpr std::array<int, 6>  ser_ih = {1, 1, 0, 2, 0, 1};
constexpr std::array<int, 6>  ser_at = {1026, 8, 0, 8, 9217, 9219};

constexpr std::array<const char*, 7>  thr = {"N", "CA", "C", "CB", "O", "CG2", "OG1"};
constexpr std::array<int, 7>  thr_ih = {1, 1, 0, 1, 0, 3, 1};
constexpr std::array<int, 7>  thr_at = {1026, 8, 0, 8, 9217, 512, 9219};

constexpr std::array<const char*, 6>  cys = {"N", "CA", "C", "CB", "O", "SG"};
constexpr std::array<int, 6>  cys_ih = {1, 1, 0, 2, 0, 1};
constexpr std::array<int, 6>  cys_at = {1026, 8, 0, 0, 9217, 9219};

constexpr std::array<const char*, 8>  met = {"N", "CA", "C", "CB", "O", "CG", "SD", "CE"};
constexpr std::array<int, 8>  met_ih = {1, 1, 0, 2, 0, 2, 0, 3};
constexpr std::array<int, 8>  met_at = {1026, 8, 0, 512, 9217, 0, 9217, 0};

constexpr std::array<const char*, 7>  pro = {"N", "CA", "C", "CB", "O", "CG", "CD"};
constexpr std::array<int, 7>  pro_ih = {0, 1, 0, 2, 0, 2, 2};
constexpr std::array<int, 7>  pro_at = {1024, 8, 0, 512, 9217, 512, 8};

constexpr std::array<const char*, 11> phe = {"N", "CA", "C", "CB", "O", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"};
constexpr std::array<int, 11> phe_ih = {1, 1, 0, 2, 0, 0, 1, 1, 1, 1, 1};
constexpr std::array<int, 11> phe_at = {1026, 8, 0, 512, 9217, 512, 512, 512, 512, 512, 512};

constexpr std::array<const char*, 12> tyr = {"N", "CA", "C", "CB", "O", "CG", "CD1", "CD2", "CE1", "CE2", "OH", "CZ"};
constexpr std::array<int, 12> tyr_ih = {1, 1, 0, 2, 0, 0, 1, 1, 1, 1, 1, 0};
constexpr std::array<int, 12> tyr_at = {1026, 8, 0, 512, 9217, 512, 512, 512, 512, 512, 9219, 0};

constexpr std::array<const char*, 14> trp = {"N", "CA", "C", "CB", "O", "CG", "CD1", "CD2", "CE2", "CE3", "NE1", "CH2", "CZ2", "CZ3"};
constexpr std::array<int, 14> trp_ih = {1, 1, 0, 2, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1};
constexpr std::array<int, 14> trp_at = {1026, 8, 0, 512, 9217, 512, 8, 512, 0, 512, 1026, 512, 512, 512};

constexpr std::array<const char*, 10> his = {"N", "CA", "C", "CB", "O", "CG", "CD2", "ND1", "CE1", "NE2"};
constexpr std::array<int, 10> his_ih = {1, 1, 0, 2, 0, 0, 1, 0, 1, 1};
constexpr std::array<int, 10> his_at = {1026, 8, 0, 512, 9217, 0, 8, 9219, 8, 9219};

constexpr std::array<const char*, 9>  glu = {"N", "CA", "C", "CB", "O", "CG", "CD", "OE1", "OE2"};
constexpr std::array<int, 9>  glu_ih = {1, 1, 0, 2, 0, 2, 0, 0, 0};
constexpr std::array<int, 9>  glu_at = {1026, 8, 0, 512, 9217, 512, 0, 9217, 9217};

constexpr std::array<const char*, 8>  asp = {"N", "CA", "C", "CB", "O", "CG", "OD1", "OD2"};
constexpr std::array<int, 8>  asp_ih = {1, 1, 0, 2, 0, 0, 0, 0};
constexpr std::array<int, 8>  asp_at = {1026, 8, 0, 512, 9217, 0, 9217, 9217};

constexpr std::array<const char*, 8>  asn = {"N", "CA", "C", "CB", "O", "CG", "ND2", "OD1"};
constexpr std::array<int, 8>  asn_ih = {1, 1, 0, 2, 0, 0, 2, 0};
constexpr std::array<int, 8>  asn_at = {1026, 8, 0, 512, 9217, 0, 1026, 9217};

constexpr std::array<const char*, 9>  gln = {"N", "CA", "C", "CB", "O", "CG", "CD", "NE2", "OE1"};
constexpr std::array<int, 9>  gln_ih = {1, 1, 0, 2, 0, 2, 0, 2, 0};
constexpr std::array<int, 9>  gln_at = {1026, 8, 0, 512, 9217, 512, 0, 1026, 9217};

constexpr std::array<const char*, 9>  lys = {"N", "CA", "C", "CB", "O", "CG", "CD", "CE", "NZ"};
constexpr std::array<int, 9>  lys_ih = {1, 1, 0, 2, 0, 2, 2, 2, 3};
constexpr std::array<int, 9>  lys_at = {1026, 8, 0, 512, 9217, 512, 512, 8, 1026};

constexpr std::array<const char*, 11> arg = {"N", "CA", "C", "CB", "O", "CG", "CD", "NE", "NH1", "NH2", "CZ"};
constexpr std::array<int, 11> arg_ih = {1, 1, 0, 2, 0, 2, 2, 1, 2, 2, 0};
constexpr std::array<int, 11> arg_at = {1026, 8, 0, 512, 9217, 512, 8, 1026, 1026, 1026, 0};

constexpr std::array<const int, 6> tyr_arom_indices = {5, 7, 9, 11, 8, 6};
constexpr std::array<const int, 6> phe_arom_indices = {5, 7, 9, 10, 8, 6};
constexpr std::array<const int, 5> his_arom_indices = {5, 7, 8, 9, 6};
constexpr std::array<const int, 5> trp_arom_indices5 = {5, 7, 8, 10, 6};
constexpr std::array<const int, 6> trp_arom_indices6 = {7, 8, 12, 11, 13, 9};


// Edge structure with canonical ordering (i < j)
struct Edge {
    int i, j;
    BondType order;

    constexpr Edge(int a, int b, BondType bt) : i(std::min(a, b)), j(std::max(a, b)), order(bt) {}

    // for sorting and binary search
    bool operator<(const Edge& other) const {
        return (i < other.i) || (i == other.i && j < other.j);
    }
};

constexpr std::array<Edge, 3> gly_edges = {
    Edge(0, 1, BondType::SINGLE),
    Edge(1, 2, BondType::SINGLE),
    Edge(2, 3, BondType::DOUBLE)
};

constexpr std::array<Edge, 4> ala_edges = {
    Edge(0, 1, BondType::SINGLE),
    Edge(1, 2, BondType::SINGLE),
    Edge(1, 3, BondType::SINGLE),
    Edge(2, 4, BondType::DOUBLE)
};

constexpr std::array<Edge, 6> val_edges = {
    Edge(0, 1, BondType::SINGLE),
    Edge(1, 2, BondType::SINGLE),
    Edge(1, 3, BondType::SINGLE),
    Edge(2, 4, BondType::DOUBLE),
    Edge(3, 5, BondType::SINGLE),
    Edge(3, 6, BondType::SINGLE)
};

constexpr std::array<Edge, 7> leu_edges = {
    Edge(0, 1, BondType::SINGLE),
    Edge(1, 2, BondType::SINGLE),
    Edge(1, 3, BondType::SINGLE),
    Edge(2, 4, BondType::DOUBLE),
    Edge(3, 5, BondType::SINGLE),
    Edge(5, 6, BondType::SINGLE),
    Edge(5, 7, BondType::SINGLE)
};

constexpr std::array<Edge, 7> ile_edges = {
    Edge(0, 1, BondType::SINGLE),
    Edge(1, 2, BondType::SINGLE),
    Edge(1, 3, BondType::SINGLE),
    Edge(2, 4, BondType::DOUBLE),
    Edge(3, 5, BondType::SINGLE),
    Edge(3, 6, BondType::SINGLE),
    Edge(5, 7, BondType::SINGLE)
};

constexpr std::array<Edge, 5> ser_edges = {
    Edge(0, 1, BondType::SINGLE),
    Edge(1, 2, BondType::SINGLE),
    Edge(1, 3, BondType::SINGLE),
    Edge(2, 4, BondType::DOUBLE),
    Edge(3, 5, BondType::SINGLE)
};

constexpr std::array<Edge, 6> thr_edges = {
    Edge(0, 1, BondType::SINGLE),
    Edge(1, 2, BondType::SINGLE),
    Edge(1, 3, BondType::SINGLE),
    Edge(2, 4, BondType::DOUBLE),
    Edge(3, 5, BondType::SINGLE),
    Edge(3, 6, BondType::SINGLE)
};

constexpr std::array<Edge, 5> cys_edges = {
    Edge(0, 1, BondType::SINGLE),
    Edge(1, 2, BondType::SINGLE),
    Edge(1, 3, BondType::SINGLE),
    Edge(2, 4, BondType::DOUBLE),
    Edge(3, 5, BondType::SINGLE)
};

constexpr std::array<Edge, 7> met_edges = {
    Edge(0, 1, BondType::SINGLE),
    Edge(1, 2, BondType::SINGLE),
    Edge(1, 3, BondType::SINGLE),
    Edge(2, 4, BondType::DOUBLE),
    Edge(3, 5, BondType::SINGLE),
    Edge(5, 6, BondType::SINGLE),
    Edge(6, 7, BondType::SINGLE)
};

constexpr std::array<Edge, 7> pro_edges = {
    Edge(0, 1, BondType::SINGLE),
    Edge(0, 6, BondType::SINGLE),
    Edge(1, 2, BondType::SINGLE),
    Edge(1, 3, BondType::SINGLE),
    Edge(2, 4, BondType::DOUBLE),
    Edge(3, 5, BondType::SINGLE),
    Edge(5, 6, BondType::SINGLE)
};

constexpr std::array<Edge, 11> phe_edges = {
    Edge(0, 1,  BondType::SINGLE),
    Edge(1, 2,  BondType::SINGLE),
    Edge(1, 3,  BondType::SINGLE),
    Edge(2, 4,  BondType::DOUBLE),
    Edge(3, 5,  BondType::SINGLE),
    Edge(5, 6,  BondType::AROMATIC),
    Edge(5, 7,  BondType::AROMATIC),
    Edge(6, 8,  BondType::AROMATIC),
    Edge(7, 9,  BondType::AROMATIC),
    Edge(8, 10, BondType::AROMATIC),
    Edge(9, 10, BondType::AROMATIC)
};

constexpr std::array<Edge, 12> tyr_edges = {
    Edge(0, 1,   BondType::SINGLE),
    Edge(1, 2,   BondType::SINGLE),
    Edge(1, 3,   BondType::SINGLE),
    Edge(2, 4,   BondType::DOUBLE),
    Edge(3, 5,   BondType::SINGLE),
    Edge(5, 6,   BondType::AROMATIC),
    Edge(5, 7,   BondType::AROMATIC),
    Edge(6, 8,   BondType::AROMATIC),
    Edge(7, 9,   BondType::AROMATIC),
    Edge(8, 11,  BondType::AROMATIC),
    Edge(9, 11,  BondType::AROMATIC),
    Edge(11, 10, BondType::SINGLE)
};

constexpr std::array<Edge, 15> trp_edges = {
    Edge(0, 1,   BondType::SINGLE),
    Edge(1, 2,   BondType::SINGLE),
    Edge(1, 3,   BondType::SINGLE),
    Edge(2, 4,   BondType::DOUBLE),
    Edge(3, 5,   BondType::SINGLE),
    Edge(5, 6,   BondType::AROMATIC),
    Edge(5, 7,   BondType::AROMATIC),
    Edge(6, 10,  BondType::AROMATIC),
    Edge(7, 8,   BondType::AROMATIC),
    Edge(7, 9,   BondType::AROMATIC),
    Edge(8, 10,  BondType::AROMATIC),
    Edge(8, 12,  BondType::AROMATIC),
    Edge(9, 13,  BondType::AROMATIC),
    Edge(12, 11, BondType::AROMATIC),
    Edge(13, 11, BondType::AROMATIC),
};

constexpr std::array<Edge, 10> his_edges = {
    Edge(0, 1, BondType::SINGLE),
    Edge(1, 2, BondType::SINGLE),
    Edge(1, 3, BondType::SINGLE),
    Edge(2, 4, BondType::DOUBLE),
    Edge(3, 5, BondType::SINGLE),
    Edge(5, 6, BondType::AROMATIC),
    Edge(5, 7, BondType::AROMATIC),
    Edge(6, 9, BondType::AROMATIC),
    Edge(7, 8, BondType::AROMATIC),
    Edge(8, 9, BondType::AROMATIC)
};

constexpr std::array<Edge, 8> glu_edges = {
    Edge(0, 1, BondType::SINGLE),
    Edge(1, 2, BondType::SINGLE),
    Edge(1, 3, BondType::SINGLE),
    Edge(2, 4, BondType::DOUBLE),
    Edge(3, 5, BondType::SINGLE),
    Edge(5, 6, BondType::SINGLE),
    Edge(6, 7, BondType::DOUBLE),
    Edge(6, 8, BondType::SINGLE)
};

constexpr std::array<Edge, 7> asp_edges = {
    Edge(0, 1, BondType::SINGLE),
    Edge(1, 2, BondType::SINGLE),
    Edge(1, 3, BondType::SINGLE),
    Edge(2, 4, BondType::DOUBLE),
    Edge(3, 5, BondType::SINGLE),
    Edge(5, 6, BondType::DOUBLE),
    Edge(5, 7, BondType::SINGLE)
};

constexpr std::array<Edge, 7> asn_edges = {
    Edge(0, 1, BondType::SINGLE),
    Edge(1, 2, BondType::SINGLE),
    Edge(1, 3, BondType::SINGLE),
    Edge(2, 4, BondType::DOUBLE),
    Edge(3, 5, BondType::SINGLE),
    Edge(5, 6, BondType::SINGLE),
    Edge(5, 7, BondType::DOUBLE)
};

constexpr std::array<Edge, 8> gln_edges = {
    Edge(0, 1, BondType::SINGLE),
    Edge(1, 2, BondType::SINGLE),
    Edge(1, 3, BondType::SINGLE),
    Edge(2, 4, BondType::DOUBLE),
    Edge(3, 5, BondType::SINGLE),
    Edge(5, 6, BondType::SINGLE),
    Edge(6, 7, BondType::SINGLE),
    Edge(6, 8, BondType::DOUBLE)
};

constexpr std::array<Edge, 8> lys_edges = {
    Edge(0, 1, BondType::SINGLE),
    Edge(1, 2, BondType::SINGLE),
    Edge(1, 3, BondType::SINGLE),
    Edge(2, 4, BondType::DOUBLE),
    Edge(3, 5, BondType::SINGLE),
    Edge(5, 6, BondType::SINGLE),
    Edge(6, 7, BondType::SINGLE),
    Edge(7, 8, BondType::SINGLE)
};

constexpr std::array<Edge, 10> arg_edges = {
    Edge(0, 1,  BondType::SINGLE),
    Edge(1, 2,  BondType::SINGLE),
    Edge(1, 3,  BondType::SINGLE),
    Edge(2, 4,  BondType::DOUBLE),
    Edge(3, 5,  BondType::SINGLE),
    Edge(5, 6,  BondType::SINGLE),
    Edge(6, 7,  BondType::SINGLE),
    Edge(7, 10, BondType::SINGLE),
    Edge(10, 8, BondType::SINGLE),
    Edge(10, 9, BondType::DOUBLE)
};

} // namespace lahuta

#endif // LAHUTA_MODEL_CONSTANTS_HPP
