#ifndef LAHUTA_BOND_TABLE_HPP
#define LAHUTA_BOND_TABLE_HPP

#include <array>
#include <gemmi/elem.hpp>

namespace lahuta {

inline constexpr int MAX_ELEMENTS = 120;
inline constexpr double dbr = 2.001f; // default bond radius

namespace types {
using AtomicNumber = int;
}

inline float bond_radius(gemmi::El el) {
  using El = gemmi::El;
  static constexpr float bond_tbl[] = {
      /*X*/ 1.42f, /*H*/ 1.42f, /*He*/ dbr,
      /*Li*/ 2.7f, /*Be*/ 2.7f, /*B*/ dbr,    /*C*/ 1.75f, /*N*/ 1.6f,   /*O*/ 1.52f,
      /*F*/ dbr,   /*Ne*/ dbr,
      /*Na*/ 2.7f, /*Mg*/ 2.7f, /*Al*/ 2.7f,  /*Si*/ 1.9f, /*P*/ 2.0f,   /*S*/ 1.9f,
      /*Cl*/ 1.8f, /*Ar*/ dbr,
      /*K*/ 2.7f,  /*Ca*/ 2.7f, /*Sc*/ 2.7f,  /*Ti*/ 2.7f, /*V*/ 2.7f,   /*Cr*/ 2.7f,
      /*Mn*/ 2.7f, /*Fe*/ 2.7f, /*Co*/ 2.7f,  /*Ni*/ 2.7f, /*Cu*/ 2.7f,  /*Zn*/ 2.7f,
      /*Ga*/ 2.7f, /*Ge*/ dbr,  /*As*/ 2.68f, /*Se*/ dbr,  /*Br*/ dbr,   /*Kr*/ dbr,
      /*Rb*/ 2.7f, /*Sr*/ 2.7f, /*Y*/ 2.7f,   /*Zr*/ 2.7f, /*Nb*/ 2.7f,  /*Mo*/ 2.7f,
      /*Tc*/ 2.7f, /*Ru*/ 2.7f, /*Rh*/ 2.7f,  /*Pd*/ 2.7f, /*Ag*/ 2.7f,  /*Cd*/ 2.7f,
      /*In*/ 2.7f, /*Sn*/ 2.7f, /*Sb*/ dbr,   /*Te*/ dbr,  /*I*/ dbr,    /*Xe*/ dbr,
      /*Cs*/ 2.7f, /*Ba*/ 2.7f, /*La*/ 2.7f,  /*Ce*/ 2.7f, /*Pr*/ 2.7f,  /*Nd*/ 2.7f,
      /*Pm*/ 2.7f, /*Sm*/ 2.7f, /*Eu*/ 2.7f,  /*Gd*/ 2.7f, /*Tb*/ 2.7f,  /*Dy*/ 2.7f,
      /*Ho*/ 2.7f, /*Er*/ 2.7f, /*Tm*/ 2.7f,  /*Yb*/ 2.7f, /*Lu*/ 2.7f,  /*Hf*/ 2.7f,
      /*Ta*/ 2.7f, /*W*/ 2.7f,  /*Re*/ 2.7f,  /*Os*/ 2.7f, /*Ir*/ 2.7f,  /*Pt*/ 2.7f,
      /*Au*/ 2.7f, /*Hg*/ 2.7f, /*Tl*/ 2.7f,  /*Pb*/ 2.7f, /*Bi*/ 2.7f,  /*Po*/ dbr,
      /*At*/ dbr,  /*Rn*/ dbr,
      /*Fr*/ 2.7f, /*Ra*/ 2.7f, /*Ac*/ 2.7f,  /*Th*/ 2.7f, /*Pa*/ 2.7f,  /*U*/ 2.7f,
      /*Np*/ 2.7f, /*Pu*/ 2.7f, /*Am*/ 2.7f,  /*Cm*/ 2.7f, /*Bk*/ 2.7f,  /*Cf*/ 2.7f,
      /*Es*/ 2.7f, /*Fm*/ 2.7f, /*Md*/ 2.7f,  /*No*/ 2.7f, /*Lr*/ 2.7f,  /*Rf*/ 2.7f,
      /*Db*/ 2.7f, /*Sg*/ 2.7f, /*Bh*/ 2.7f,  /*Hs*/ 2.7f, /*Mt*/ 2.88f, /*Ds*/ dbr,
      /*Rg*/ dbr,  /*Cn*/ dbr,  /*Nh*/ dbr,   /*Fl*/ dbr,  /*Mc*/ dbr,   /*Lv*/ dbr,
      /*Ts*/ dbr,  /*Og*/ dbr,
      /*D*/ 1.42f, /*END*/ 0.0f};

  static_assert(bond_tbl[static_cast<int>(El::D)] == 1.42f, "Error in bond radius for Deuterium");
  static_assert(sizeof(bond_tbl) / sizeof(bond_tbl[0]) == static_cast<int>(El::END) + 1, "Size mismatch");
  return bond_tbl[static_cast<int>(el)];
}

inline constexpr std::array<std::array<double, MAX_ELEMENTS>, MAX_ELEMENTS> PairThresholds = []() {
  constexpr std::array<std::tuple<int, int, double>, 63> pair_values{
      {{1, 1, 0.8},    {1, 5, 1.31},   {1, 6, 1.2},   {1, 7, 1.15},   {1, 8, 1.1},    {1, 9, 1.0},
       {5, 5, 1.84},   {5, 6, 1.88},   {6, 6, 1.75},  {5, 7, 1.56},   {4, 8, 1.76},   {6, 7, 1.6},
       {5, 8, 1.68},   {4, 9, 1.63},   {7, 7, 1.6},   {6, 8, 1.59},   {5, 9, 1.36},   {6, 9, 1.45},
       {1, 15, 1.47},  {8, 8, 1.6},    {1, 16, 1.45}, {1, 17, 1.4},   {9, 9, 1.55},   {7, 12, 2.4},
       {8, 12, 2.24},  {6, 14, 1.91},  {5, 15, 1.98}, {9, 12, 2.02},  {6, 16, 2.0},   {6, 17, 1.9},
       {8, 16, 1.8},   {14, 14, 2.37}, {15, 15, 2.3}, {15, 16, 2.3},  {16, 16, 2.3},  {17, 17, 2.1},
       {1, 34, 1.54},  {1, 35, 1.0},   {6, 33, 2.6},  {6, 34, 2.27},  {8, 33, 1.93},  {6, 35, 2.1},
       {8, 34, 2.05},  {7, 35, 2.06},  {8, 35, 1.62}, {16, 33, 2.68}, {16, 34, 2.33}, {1, 53, 1.0},
       {6, 52, 2.14},  {6, 53, 2.48},  {8, 52, 2.1},  {8, 53, 1.72},  {34, 34, 2.34}, {35, 46, 2.44},
       {7, 78, 2.11},  {8, 78, 2.6},   {6, 80, 2.36}, {16, 80, 2.75}, {53, 53, 2.73}, {35, 73, 2.63},
       {35, 78, 2.84}, {35, 80, 2.87}, {53, 80, 2.81}}};

  std::array<std::array<double, MAX_ELEMENTS>, MAX_ELEMENTS> thresholds{};
  for (const auto &[i, j, value] : pair_values) {
    thresholds[i][j] = value;
  }

  // mirror values across the diagonal
  for (int i = 1; i < MAX_ELEMENTS; ++i) {
    for (int j = i + 1; j < MAX_ELEMENTS; ++j) {
      thresholds[j][i] = thresholds[i][j];
    }
  }

  return thresholds;
}();

/// retrieve the bond radius for an element
inline double get_bond_radius(types::AtomicNumber i) {
  gemmi::El e = static_cast<gemmi::El>(i);
  if (e < gemmi::El::H || e >= gemmi::El::END) return dbr;
  return bond_radius(e);
}

/// retrieve the precomputed pair threshold for two elements
inline constexpr double get_element_pair_threshold(types::AtomicNumber i, types::AtomicNumber j) {
  // RDKit stores atomic numbers as unsigned integers (std::uint8_t)
  if (i >= MAX_ELEMENTS || j >= MAX_ELEMENTS) return 0.0;
  return PairThresholds[i][j];
}

/// retrieve the pair threshold for two elements
inline double get_pair_threshold(types::AtomicNumber a, types::AtomicNumber b) {
  double threshold_ab = get_element_pair_threshold(a, b);
  if (threshold_ab > 0.0) return threshold_ab;

  return (get_bond_radius(a) + get_bond_radius(b)) / 1.95;
}

} // namespace lahuta

#endif // LAHUTA_BOND_TABLE_HPP
