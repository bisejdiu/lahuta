/* C++ code produced by gperf version 3.1 */
/* Command-line: gperf --constants-prefix=T_ -L C++ -t -N getToken -K token -H
 * hashToken -Z TokenTable tokens.gperf  */
/* Computed positions: -k'1-3' */

#if !(                                                                         \
    (' ' == 32) && ('!' == 33) && ('"' == 34) && ('#' == 35) && ('%' == 37) && \
    ('&' == 38) && ('\'' == 39) && ('(' == 40) && (')' == 41) &&               \
    ('*' == 42) && ('+' == 43) && (',' == 44) && ('-' == 45) && ('.' == 46) && \
    ('/' == 47) && ('0' == 48) && ('1' == 49) && ('2' == 50) && ('3' == 51) && \
    ('4' == 52) && ('5' == 53) && ('6' == 54) && ('7' == 55) && ('8' == 56) && \
    ('9' == 57) && (':' == 58) && (';' == 59) && ('<' == 60) && ('=' == 61) && \
    ('>' == 62) && ('?' == 63) && ('A' == 65) && ('B' == 66) && ('C' == 67) && \
    ('D' == 68) && ('E' == 69) && ('F' == 70) && ('G' == 71) && ('H' == 72) && \
    ('I' == 73) && ('J' == 74) && ('K' == 75) && ('L' == 76) && ('M' == 77) && \
    ('N' == 78) && ('O' == 79) && ('P' == 80) && ('Q' == 81) && ('R' == 82) && \
    ('S' == 83) && ('T' == 84) && ('U' == 85) && ('V' == 86) && ('W' == 87) && \
    ('X' == 88) && ('Y' == 89) && ('Z' == 90) && ('[' == 91) &&                \
    ('\\' == 92) && (']' == 93) && ('^' == 94) && ('_' == 95) &&               \
    ('a' == 97) && ('b' == 98) && ('c' == 99) && ('d' == 100) &&               \
    ('e' == 101) && ('f' == 102) && ('g' == 103) && ('h' == 104) &&            \
    ('i' == 105) && ('j' == 106) && ('k' == 107) && ('l' == 108) &&            \
    ('m' == 109) && ('n' == 110) && ('o' == 111) && ('p' == 112) &&            \
    ('q' == 113) && ('r' == 114) && ('s' == 115) && ('t' == 116) &&            \
    ('u' == 117) && ('v' == 118) && ('w' == 119) && ('x' == 120) &&            \
    ('y' == 121) && ('z' == 122) && ('{' == 123) && ('|' == 124) &&            \
    ('}' == 125) && ('~' == 126))
/* The character set is not based on ISO-646.  */
#error                                                                         \
    "gperf generated tables don't work with this execution character set. Please report a bug to <bug-gperf@gnu.org>."
#endif

#include <cstring>
struct TokenNames {
  const char *token;
};

#define T_TOTAL_KEYWORDS 61
#define T_MIN_WORD_LENGTH 1
#define T_MAX_WORD_LENGTH 4
#define T_MIN_HASH_VALUE 1
#define T_MAX_HASH_VALUE 128
/* maximum key range = 128, duplicates = 0 */

class TokenTable {
private:
  static inline unsigned int hashToken(const char *str, size_t len);

public:
  static struct TokenNames *getToken(const char *str, size_t len);
};

inline unsigned int TokenTable::hashToken(const char *str, size_t len) {
  static unsigned char asso_values[] = {
      129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129,
      129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129,
      129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129,
      129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129,
      129, 129, 129, 129, 129, 0,   2,   40,  0,   25,  26,  10,  30,  15,  129,
      129, 25,  11,  55,  1,   45,  10,  7,   5,   30,  5,   56,  129, 129, 55,
      16,  129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129,
      129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129,
      129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129,
      129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129,
      129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129,
      129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129,
      129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129,
      129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129,
      129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129,
      129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129,
      129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129,
      129, 129};
  unsigned int hval = len;

  switch (hval) {
  default:
    hval += asso_values[static_cast<unsigned char>(str[2] + 1)];
  /*FALLTHROUGH*/
  case 2:
    hval += asso_values[static_cast<unsigned char>(str[1])];
  /*FALLTHROUGH*/
  case 1:
    hval += asso_values[static_cast<unsigned char>(str[0])];
    break;
  }
  return hval;
}

inline struct TokenNames *TokenTable::getToken(const char *str, size_t len) {
  static struct TokenNames wordlist[] = {
      {""},

      {"A"},

      {"DA"},   {""}, {""}, {""},

      {"U"},

      {"DU"},

      {"DAR"},

      {"ASN"},  {""},

      {"G"},

      {"DG"},   {""},

      {"DAL"},  {""},

      {"I"},

      {"DI"},

      {"ASP"},

      {"ASPP"}, {""}, {""}, {""},

      {"ASH"},  {""}, {""}, {""}, {""}, {""}, {""},

      {"ALA"},

      {"T"},

      {"DT"},

      {"SEC"},

      {"LSN"},  {""}, {""}, {""},

      {"SER"},

      {"GLN"},

      {"ARG"},

      {"C"},

      {"DC"},

      {"SEP"},

      {"MET"},

      {"MSE"},  {""}, {""},

      {"HSP"},

      {"APN"},

      {"TRP"},  {""}, {""},

      {"GLH"},

      {"GLY"},  {""},

      {"N"},

      {"DN"},

      {"HIP"},

      {"GPN"},  {""}, {""}, {""},

      {"HSD"},

      {"HSE"},  {""}, {""}, {""},

      {"THR"},

      {"ILE"},

      {"VAL"},  {""}, {""},

      {"HID"},

      {"HIE"},  {""}, {""}, {""},

      {"HIS"},

      {"TPN"},  {""}, {""}, {""},

      {"PTR"},

      {"LYN"},  {""}, {""}, {""},

      {"UNK"},

      {"CPN"},

      {"PCA"},  {""}, {""},

      {"TYR"},

      {"GLU"},

      {"GLUP"}, {""}, {""},

      {"HYP"},  {""},

      {"PRO"},  {""}, {""}, {""},

      {"PHE"},  {""}, {""}, {""}, {""},

      {"LEU"},  {""}, {""}, {""},

      {"LYS"},

      {"PYL"},  {""}, {""}, {""}, {""}, {""}, {""}, {""}, {""},

      {"TPO"},  {""}, {""}, {""}, {""},

      {"CYS"}};

  if (len <= T_MAX_WORD_LENGTH && len >= T_MIN_WORD_LENGTH) {
    unsigned int key = hashToken(str, len);

    if (key <= T_MAX_HASH_VALUE) {
      const char *s = wordlist[key].token;

      if (*str == *s && !strcmp(str + 1, s + 1))
        return &wordlist[key];
    }
  }
  return 0;
}

