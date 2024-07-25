/* C++ code produced by gperf version 3.1 */
/* Command-line: gperf --constants-prefix=RN_ -L C++ -t -N getResName -K resName -H hashResName -Z ResNameTable aa_names.gperf  */
/* Computed positions: -k'1-3' */

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

#line 1 "aa_names.gperf"

#include <cstring>
#line 5 "aa_names.gperf"
struct ProteinResNames {
    const char* resName;
};

#define RN_TOTAL_KEYWORDS 43
#define RN_MIN_WORD_LENGTH 3
#define RN_MAX_WORD_LENGTH 4
#define RN_MIN_HASH_VALUE 3
#define RN_MAX_HASH_VALUE 154
/* maximum key range = 152, duplicates = 0 */

class ResNameTable
{
private:
  static inline unsigned int hashResName (const char *str, size_t len);
public:
  static struct ProteinResNames *getResName (const char *str, size_t len);
};

inline unsigned int
ResNameTable::hashResName (const char *str, size_t len)
{
  static unsigned char asso_values[] =
    {
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155,  30,  55,  30,   1,   0,
       45,  50,  15,   5, 155, 155,  40,  11,  11,   6,
       25,  10,  21,   0,  30,   0,  60, 155, 155,  45,
       30, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155, 155, 155, 155,
      155, 155, 155, 155, 155, 155, 155
    };
  return len + asso_values[static_cast<unsigned char>(str[2]+1)] + asso_values[static_cast<unsigned char>(str[1])] + asso_values[static_cast<unsigned char>(str[0])];
}

inline struct ProteinResNames *
ResNameTable::getResName (const char *str, size_t len)
{
  static struct ProteinResNames wordlist[] =
    {
      {""}, {""}, {""},
#line 24 "aa_names.gperf"
      {"SER"},
#line 30 "aa_names.gperf"
      {"SEC"},
      {""}, {""}, {""}, {""}, {""}, {""}, {""}, {""},
#line 34 "aa_names.gperf"
      {"SEP"},
#line 18 "aa_names.gperf"
      {"MET"},
      {""}, {""}, {""},
#line 39 "aa_names.gperf"
      {"HSD"},
      {""}, {""}, {""}, {""},
#line 45 "aa_names.gperf"
      {"HID"},
      {""}, {""}, {""}, {""},
#line 41 "aa_names.gperf"
      {"HSP"},
      {""}, {""}, {""}, {""},
#line 47 "aa_names.gperf"
      {"HIP"},
#line 52 "aa_names.gperf"
      {"DAR"},
      {""}, {""}, {""},
#line 49 "aa_names.gperf"
      {"ASH"},
#line 21 "aa_names.gperf"
      {"ASN"},
      {""}, {""}, {""},
#line 27 "aa_names.gperf"
      {"ASP"},
#line 43 "aa_names.gperf"
      {"ASPP"},
#line 51 "aa_names.gperf"
      {"DAL"},
      {""}, {""},
#line 29 "aa_names.gperf"
      {"THR"},
#line 42 "aa_names.gperf"
      {"LSN"},
      {""}, {""}, {""},
#line 10 "aa_names.gperf"
      {"HIS"},
#line 32 "aa_names.gperf"
      {"UNK"},
      {""}, {""}, {""},
#line 36 "aa_names.gperf"
      {"PTR"},
#line 33 "aa_names.gperf"
      {"MSE"},
      {""}, {""}, {""},
#line 40 "aa_names.gperf"
      {"HSE"},
#line 16 "aa_names.gperf"
      {"TRP"},
      {""}, {""}, {""},
#line 46 "aa_names.gperf"
      {"HIE"},
#line 11 "aa_names.gperf"
      {"ARG"},
      {""}, {""}, {""},
#line 38 "aa_names.gperf"
      {"HYP"},
#line 19 "aa_names.gperf"
      {"PRO"},
      {""}, {""}, {""},
#line 26 "aa_names.gperf"
      {"TYR"},
      {""}, {""}, {""}, {""},
#line 35 "aa_names.gperf"
      {"TPO"},
#line 31 "aa_names.gperf"
      {"PYL"},
      {""}, {""}, {""},
#line 14 "aa_names.gperf"
      {"PHE"},
      {""}, {""}, {""}, {""},
#line 13 "aa_names.gperf"
      {"ILE"},
#line 48 "aa_names.gperf"
      {"LYN"},
      {""}, {""}, {""},
#line 50 "aa_names.gperf"
      {"GLH"},
#line 25 "aa_names.gperf"
      {"GLN"},
      {""}, {""}, {""},
#line 15 "aa_names.gperf"
      {"LEU"},
#line 22 "aa_names.gperf"
      {"VAL"},
      {""}, {""}, {""},
#line 20 "aa_names.gperf"
      {"CYS"},
      {""}, {""}, {""}, {""},
#line 37 "aa_names.gperf"
      {"PCA"},
      {""}, {""}, {""}, {""},
#line 12 "aa_names.gperf"
      {"LYS"},
      {""}, {""}, {""}, {""},
#line 23 "aa_names.gperf"
      {"GLY"},
      {""}, {""}, {""}, {""},
#line 17 "aa_names.gperf"
      {"ALA"},
      {""}, {""}, {""}, {""}, {""}, {""}, {""}, {""}, {""},
      {""}, {""}, {""}, {""}, {""}, {""}, {""}, {""}, {""},
      {""}, {""}, {""}, {""}, {""}, {""},
#line 28 "aa_names.gperf"
      {"GLU"},
#line 44 "aa_names.gperf"
      {"GLUP"}
    };

  if (len <= RN_MAX_WORD_LENGTH && len >= RN_MIN_WORD_LENGTH)
    {
      unsigned int key = hashResName (str, len);

      if (key <= RN_MAX_HASH_VALUE)
        {
          const char *s = wordlist[key].resName;

          if (*str == *s && !strcmp (str + 1, s + 1))
            return &wordlist[key];
        }
    }
  return 0;
}
#line 53 "aa_names.gperf"

