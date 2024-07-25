/* C++ code produced by gperf version 3.1 */
/* Command-line: gperf -L C++ -t -N get_bond_order -K order -H hash_func -Z BondOrderTable bond_order3.gperf  */
/* Computed positions: -k'1-3,6,$' */

#if !((' ' == 32) && ('!' == 33) && ('"' == 34) && ('#' == 35) \
      && ('%' == 37) && ('&' == 38) && ('\'' == 39) && ('(' == 40) \
      && (')' == 41) && ('*' == 42) && ('+' == 43) && (',' == 44) \
      && ('-' == 45) && ('.' == 46) && ('/' == 47) && ('0' == 48) \
      && ('1' == 49) && ('2' == 50) && ('3' == 51) && ('4' == 52) \
      && ('5' == 53) && ('6' == 54) && ('7' == 55) && ('8' == 56) \
      && ('9' == 57) && (':' == 58) && (';' == 59) && ('<' == 60) \
      && ('=' == 61) && ('>' == 62) && ('?' == 63) && ('A' == 65) \
      && ('B' == 66) && ('C' == 67) && ('D' == 68) && ('E' == 69) \
      && ('F' == 70) && ('G' == 71) && ('H' == 72) && ('I' == 73) \
      && ('J' == 74) && ('K' == 75) && ('L' == 76) && ('M' == 77) \
      && ('N' == 78) && ('O' == 79) && ('P' == 80) && ('Q' == 81) \
      && ('R' == 82) && ('S' == 83) && ('T' == 84) && ('U' == 85) \
      && ('V' == 86) && ('W' == 87) && ('X' == 88) && ('Y' == 89) \
      && ('Z' == 90) && ('[' == 91) && ('\\' == 92) && (']' == 93) \
      && ('^' == 94) && ('_' == 95) && ('a' == 97) && ('b' == 98) \
      && ('c' == 99) && ('d' == 100) && ('e' == 101) && ('f' == 102) \
      && ('g' == 103) && ('h' == 104) && ('i' == 105) && ('j' == 106) \
      && ('k' == 107) && ('l' == 108) && ('m' == 109) && ('n' == 110) \
      && ('o' == 111) && ('p' == 112) && ('q' == 113) && ('r' == 114) \
      && ('s' == 115) && ('t' == 116) && ('u' == 117) && ('v' == 118) \
      && ('w' == 119) && ('x' == 120) && ('y' == 121) && ('z' == 122) \
      && ('{' == 123) && ('|' == 124) && ('}' == 125) && ('~' == 126))
/* The character set is not based on ISO-646.  */
#error "gperf generated tables don't work with this execution character set. Please report a bug to <bug-gperf@gnu.org>."
#endif

#line 1 "bond_order3.gperf"

#include <cstring>
#line 5 "bond_order3.gperf"
struct BondOrder {
    const char* order;
};

#define TOTAL_KEYWORDS 45
#define MIN_WORD_LENGTH 7
#define MAX_WORD_LENGTH 11
#define MIN_HASH_VALUE 8
#define MAX_HASH_VALUE 111
/* maximum key range = 104, duplicates = 0 */

class BondOrderTable
{
private:
  static inline unsigned int hash_func (const char *str, size_t len);
public:
  static struct BondOrder *get_bond_order (const char *str, size_t len);
};

inline unsigned int
BondOrderTable::hash_func (const char *str, size_t len)
{
  static unsigned char asso_values[] =
    {
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112,  46,
       30,  10,  45,  25,  35,   0, 112, 112, 112, 112,
      112, 112, 112, 112, 112,   5, 112,  10,   0,   0,
      112,   0,   5,  10, 112, 112,  35, 112,   0,  10,
        5, 112,   0,  15,  20,  20, 112, 112, 112,  25,
        5, 112, 112, 112, 112,   0, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112, 112, 112, 112, 112,
      112, 112, 112, 112, 112, 112
    };
  return len + asso_values[static_cast<unsigned char>(str[5])] + asso_values[static_cast<unsigned char>(str[2])] + asso_values[static_cast<unsigned char>(str[1])] + asso_values[static_cast<unsigned char>(str[0])] + asso_values[static_cast<unsigned char>(str[len - 1])];
}


inline BondOrder *
BondOrderTable::get_bond_order (const char *str, size_t len)
{
  static struct BondOrder wordlist[] =
    {
      {""}, {""}, {""}, {""}, {""}, {""}, {""}, {""},
#line 41 "bond_order3.gperf"
      {"DG_C8_N7"},
      {""}, {""}, {""}, {""},
#line 51 "bond_order3.gperf"
      {"DA_C8_N7"},
      {""}, {""}, {""},
#line 27 "bond_order3.gperf"
      {"G_C8_N7"},
#line 43 "bond_order3.gperf"
      {"DG_C2_N3"},
      {""},
#line 15 "bond_order3.gperf"
      {"PHE_CD1_CG"},
      {""},
#line 37 "bond_order3.gperf"
      {"A_C8_N7"},
#line 48 "bond_order3.gperf"
      {"DA_C2_N3"},
      {""},
#line 13 "bond_order3.gperf"
      {"PHE_CE1_CZ"},
      {""},
#line 29 "bond_order3.gperf"
      {"G_C2_N3"},
#line 45 "bond_order3.gperf"
      {"DC_C4_N3"},
      {""}, {""}, {""},
#line 34 "bond_order3.gperf"
      {"A_C2_N3"},
#line 42 "bond_order3.gperf"
      {"DG_C4_C5"},
      {""},
#line 16 "bond_order3.gperf"
      {"TRP_CD1_CG"},
      {""},
#line 31 "bond_order3.gperf"
      {"C_C4_N3"},
#line 50 "bond_order3.gperf"
      {"DA_C4_C5"},
      {""},
#line 10 "bond_order3.gperf"
      {"HIS_CD2_CG"},
      {""}, {""},
#line 44 "bond_order3.gperf"
      {"DG_C6_O6"},
      {""}, {""},
#line 18 "bond_order3.gperf"
      {"TRP_CE3_CZ3"},
      {""},
#line 47 "bond_order3.gperf"
      {"DC_C2_O2"},
      {""},
#line 12 "bond_order3.gperf"
      {"ARG_CZ_NH2"},
#line 14 "bond_order3.gperf"
      {"PHE_CD2_CE2"},
#line 28 "bond_order3.gperf"
      {"G_C4_C5"},
#line 46 "bond_order3.gperf"
      {"DC_C5_C6"},
      {""},
#line 22 "bond_order3.gperf"
      {"TYR_CD1_CG"},
      {""},
#line 36 "bond_order3.gperf"
      {"A_C4_C5"},
#line 53 "bond_order3.gperf"
      {"DT_C2_O2"},
#line 49 "bond_order3.gperf"
      {"DA_C6_N1"},
#line 24 "bond_order3.gperf"
      {"TYR_CE1_CZ"},
      {""},
#line 30 "bond_order3.gperf"
      {"G_C6_O6"},
#line 52 "bond_order3.gperf"
      {"DT_C5_C6"},
      {""}, {""},
#line 17 "bond_order3.gperf"
      {"TRP_CD2_CE2"},
#line 33 "bond_order3.gperf"
      {"C_C2_O2"},
#line 35 "bond_order3.gperf"
      {"A_C6_N1"},
      {""}, {""},
#line 19 "bond_order3.gperf"
      {"TRP_CH2_CZ2"},
#line 32 "bond_order3.gperf"
      {"C_C5_C6"},
#line 54 "bond_order3.gperf"
      {"DT_C4_O4"},
      {""}, {""},
#line 20 "bond_order3.gperf"
      {"ASN_CG_OD1"},
#line 39 "bond_order3.gperf"
      {"U_C2_O2"},
      {""}, {""}, {""},
#line 25 "bond_order3.gperf"
      {"ASP_CG_OD1"},
#line 38 "bond_order3.gperf"
      {"U_C5_C6"},
      {""}, {""}, {""},
#line 23 "bond_order3.gperf"
      {"TYR_CD2_CE2"},
#line 11 "bond_order3.gperf"
      {"HIS_CE1_ND1"},
      {""}, {""}, {""},
#line 21 "bond_order3.gperf"
      {"GLN_CD_OE1"},
#line 40 "bond_order3.gperf"
      {"U_C4_O4"},
      {""}, {""}, {""}, {""}, {""}, {""}, {""}, {""}, {""},
      {""}, {""}, {""}, {""}, {""}, {""}, {""}, {""}, {""},
#line 26 "bond_order3.gperf"
      {"GLU_CD_OE1"}
    };

  if (len <= MAX_WORD_LENGTH && len >= MIN_WORD_LENGTH)
    {
      unsigned int key = hash_func (str, len);

      if (key <= MAX_HASH_VALUE)
        {
          const char *s = wordlist[key].order;

          if (*str == *s && !strcmp (str + 1, s + 1))
            return &wordlist[key];
        }
    }
  return 0;
}
#line 55 "bond_order3.gperf"

