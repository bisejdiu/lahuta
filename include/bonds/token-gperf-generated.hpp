#ifndef _TOKEN_H_
#define _TOKEN_H_

/* C++ code produced by gperf version 3.1 */
/* Command-line: gperf token.gperf  */
/* Computed positions: -k'1-2,$' */

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

#line 12 "token.gperf"

#include <stddef.h>
#include <string.h>
#include "token.h"

namespace lahuta {
namespace {
#line 20 "token.gperf"
struct keywordEntry {
  int string_offset;
  resTokenType type;
};
/* maximum key range = 333, duplicates = 0 */

class resNameKeyword
{
private:
  static inline unsigned int hash (const char *str, size_t len);
public:
  static const struct keywordEntry *look_up (const char *str, size_t len);
};

inline unsigned int
resNameKeyword::hash (const char *str, size_t len)
{
  static const unsigned short asso_values[] =
    {
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      100,  15,   5, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334,  25,  32,  20,  10,  90,
      115,  15,   5,  50,   2, 334,   7,  50,  65, 100,
       30,   2, 110, 105,   0,  90,  15,  27, 334, 117,
        7, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334, 334, 334, 334,
      334, 334, 334, 334, 334, 334, 334
    };
  unsigned int hval = len;

  switch (hval)
    {
      default:
        hval += asso_values[static_cast<unsigned char>(str[1]+1)];
      /*FALLTHROUGH*/
      case 1:
        hval += asso_values[static_cast<unsigned char>(str[0])];
        break;
    }
  return hval + asso_values[static_cast<unsigned char>(str[len - 1])];
}

struct stringpool_t
  {
    char stringpool_str1[sizeof("T")];
    char stringpool_str11[sizeof("TIP4")];
    char stringpool_str18[sizeof("HSD")];
    char stringpool_str20[sizeof("HID")];
    char stringpool_str21[sizeof("TIP3")];
    char stringpool_str22[sizeof("DIL")];
    char stringpool_str25[sizeof("DGL")];
    char stringpool_str28[sizeof("DSG")];
    char stringpool_str31[sizeof("G")];
    char stringpool_str32[sizeof("DG")];
    char stringpool_str33[sizeof("ASH")];
    char stringpool_str35[sizeof("TIP")];
    char stringpool_str38[sizeof("HSP")];
    char stringpool_str40[sizeof("HIP")];
    char stringpool_str41[sizeof("C")];
    char stringpool_str42[sizeof("DC")];
    char stringpool_str43[sizeof("HOH")];
    char stringpool_str45[sizeof("HYP")];
    char stringpool_str47[sizeof("PYL")];
    char stringpool_str51[sizeof("A")];
    char stringpool_str52[sizeof("DAL")];
    char stringpool_str53[sizeof("DOD")];
    char stringpool_str55[sizeof("W")];
    char stringpool_str57[sizeof("VAL")];
    char stringpool_str58[sizeof("ASP")];
    char stringpool_str59[sizeof("ASPP")];
    char stringpool_str62[sizeof("WAT")];
    char stringpool_str63[sizeof("HOA")];
    char stringpool_str64[sizeof("DI")];
    char stringpool_str65[sizeof("DVA")];
    char stringpool_str68[sizeof("PCA")];
    char stringpool_str69[sizeof("DA")];
    char stringpool_str70[sizeof("TPN")];
    char stringpool_str73[sizeof("GLH")];
    char stringpool_str75[sizeof("LSN")];
    char stringpool_str78[sizeof("DSN")];
    char stringpool_str80[sizeof("DPN")];
    char stringpool_str82[sizeof("LYN")];
    char stringpool_str83[sizeof("DGN")];
    char stringpool_str85[sizeof("GPN")];
    char stringpool_str88[sizeof("MOH")];
    char stringpool_str90[sizeof("CPN")];
    char stringpool_str93[sizeof("ASN")];
    char stringpool_str95[sizeof("APN")];
    char stringpool_str98[sizeof("HSE")];
    char stringpool_str99[sizeof("GLUP")];
    char stringpool_str100[sizeof("HIE")];
    char stringpool_str101[sizeof("I")];
    char stringpool_str102[sizeof("DT")];
    char stringpool_str103[sizeof("ALA")];
    char stringpool_str105[sizeof("TPO")];
    char stringpool_str108[sizeof("DTH")];
    char stringpool_str113[sizeof("DHI")];
    char stringpool_str115[sizeof("HIS")];
    char stringpool_str117[sizeof("DU")];
    char stringpool_str118[sizeof("D3O")];
    char stringpool_str120[sizeof("TYR")];
    char stringpool_str122[sizeof("LYS")];
    char stringpool_str123[sizeof("H2O")];
    char stringpool_str125[sizeof("DPR")];
    char stringpool_str128[sizeof("ACE")];
    char stringpool_str130[sizeof("SPC")];
    char stringpool_str131[sizeof("N")];
    char stringpool_str133[sizeof("GLN")];
    char stringpool_str135[sizeof("CYS")];
    char stringpool_str138[sizeof("TRP")];
    char stringpool_str140[sizeof("DCY")];
    char stringpool_str143[sizeof("MSE")];
    char stringpool_str145[sizeof("SOL")];
    char stringpool_str148[sizeof("ARG")];
    char stringpool_str150[sizeof("DAS")];
    char stringpool_str153[sizeof("DLE")];
    char stringpool_str155[sizeof("DAR")];
    char stringpool_str158[sizeof("GLU")];
    char stringpool_str163[sizeof("THR")];
    char stringpool_str168[sizeof("MET")];
    char stringpool_str173[sizeof("PHE")];
    char stringpool_str177[sizeof("DN")];
    char stringpool_str178[sizeof("MED")];
    char stringpool_str180[sizeof("DLY")];
    char stringpool_str181[sizeof("U")];
    char stringpool_str183[sizeof("FMT")];
    char stringpool_str185[sizeof("GLY")];
    char stringpool_str188[sizeof("NEH")];
    char stringpool_str193[sizeof("ILE")];
    char stringpool_str198[sizeof("E1H")];
    char stringpool_str203[sizeof("DNE")];
    char stringpool_str213[sizeof("DTR")];
    char stringpool_str215[sizeof("LEU")];
    char stringpool_str218[sizeof("NH2")];
    char stringpool_str220[sizeof("DTY")];
    char stringpool_str223[sizeof("NME")];
    char stringpool_str233[sizeof("PTR")];
    char stringpool_str238[sizeof("PRO")];
    char stringpool_str243[sizeof("SEC")];
    char stringpool_str253[sizeof("SEP")];
    char stringpool_str258[sizeof("FOR")];
    char stringpool_str333[sizeof("SER")];
  };
static const struct stringpool_t stringpool_contents =
  {
    "T",
    "TIP4",
    "HSD",
    "HID",
    "TIP3",
    "DIL",
    "DGL",
    "DSG",
    "G",
    "DG",
    "ASH",
    "TIP",
    "HSP",
    "HIP",
    "C",
    "DC",
    "HOH",
    "HYP",
    "PYL",
    "A",
    "DAL",
    "DOD",
    "W",
    "VAL",
    "ASP",
    "ASPP",
    "WAT",
    "HOA",
    "DI",
    "DVA",
    "PCA",
    "DA",
    "TPN",
    "GLH",
    "LSN",
    "DSN",
    "DPN",
    "LYN",
    "DGN",
    "GPN",
    "MOH",
    "CPN",
    "ASN",
    "APN",
    "HSE",
    "GLUP",
    "HIE",
    "I",
    "DT",
    "ALA",
    "TPO",
    "DTH",
    "DHI",
    "HIS",
    "DU",
    "D3O",
    "TYR",
    "LYS",
    "H2O",
    "DPR",
    "ACE",
    "SPC",
    "N",
    "GLN",
    "CYS",
    "TRP",
    "DCY",
    "MSE",
    "SOL",
    "ARG",
    "DAS",
    "DLE",
    "DAR",
    "GLU",
    "THR",
    "MET",
    "PHE",
    "DN",
    "MED",
    "DLY",
    "U",
    "FMT",
    "GLY",
    "NEH",
    "ILE",
    "E1H",
    "DNE",
    "DTR",
    "LEU",
    "NH2",
    "DTY",
    "NME",
    "PTR",
    "PRO",
    "SEC",
    "SEP",
    "FOR",
    "SER"
  };
#define stringpool ((const char *) &stringpool_contents)
const struct keywordEntry *
resNameKeyword::look_up (const char *str, size_t len)
{
  enum
    {
      TOTAL_KEYWORDS = 98,
      MIN_WORD_LENGTH = 1,
      MAX_WORD_LENGTH = 4,
      MIN_HASH_VALUE = 1,
      MAX_HASH_VALUE = 333
    };

  static const struct keywordEntry wordlist[] =
    {
      {-1,static_cast<resTokenType>(0)},
#line 89 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str1,      resTokenType::T},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 122 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str11,   resTokenType::TIP4},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 67 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str18,    resTokenType::HSD},
      {-1,static_cast<resTokenType>(0)},
#line 70 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str20,    resTokenType::HID},
#line 121 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str21,   resTokenType::TIP3},
#line 55 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str22,    resTokenType::DIL},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 52 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str25,    resTokenType::DGL},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 49 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str28,    resTokenType::DSG},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 88 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str31,      resTokenType::G},
#line 95 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str32,     resTokenType::DG},
#line 74 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str33,    resTokenType::ASH},
      {-1,static_cast<resTokenType>(0)},
#line 120 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str35,    resTokenType::TIP},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 69 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str38,    resTokenType::HSP},
      {-1,static_cast<resTokenType>(0)},
#line 72 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str40,    resTokenType::HIP},
#line 87 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str41,      resTokenType::C},
#line 94 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str42,     resTokenType::DC},
#line 115 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str43,    resTokenType::HOH},
      {-1,static_cast<resTokenType>(0)},
#line 83 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str45,    resTokenType::HYP},
      {-1,static_cast<resTokenType>(0)},
#line 85 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str47,    resTokenType::PYL},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 86 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str51,      resTokenType::A},
#line 47 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str52,    resTokenType::DAL},
#line 118 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str53,    resTokenType::DOD},
      {-1,static_cast<resTokenType>(0)},
#line 117 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str55,      resTokenType::W},
      {-1,static_cast<resTokenType>(0)},
#line 28 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str57,    resTokenType::VAL},
#line 41 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str58,    resTokenType::ASP},
#line 75 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str59,   resTokenType::ASPP},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 114 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str62,    resTokenType::WAT},
#line 110 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str63,    resTokenType::HOA},
#line 97 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str64,     resTokenType::DI},
#line 65 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str65,    resTokenType::DVA},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 46 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str68,    resTokenType::PCA},
#line 93 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str69,     resTokenType::DA},
#line 102 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str70,    resTokenType::TPN},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 76 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str73,    resTokenType::GLH},
      {-1,static_cast<resTokenType>(0)},
#line 77 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str75,    resTokenType::LSN},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 61 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str78,    resTokenType::DSN},
      {-1,static_cast<resTokenType>(0)},
#line 59 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str80,    resTokenType::DPN},
      {-1,static_cast<resTokenType>(0)},
#line 78 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str82,    resTokenType::LYN},
#line 53 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str83,    resTokenType::DGN},
      {-1,static_cast<resTokenType>(0)},
#line 103 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str85,    resTokenType::GPN},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 112 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str88,    resTokenType::MOH},
      {-1,static_cast<resTokenType>(0)},
#line 101 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str90,    resTokenType::CPN},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 42 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str93,    resTokenType::ASN},
      {-1,static_cast<resTokenType>(0)},
#line 100 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str95,    resTokenType::APN},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 68 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str98,    resTokenType::HSE},
#line 73 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str99,   resTokenType::GLUP},
#line 71 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str100,    resTokenType::HIE},
#line 91 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str101,      resTokenType::I},
#line 96 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str102,     resTokenType::DT},
#line 27 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str103,    resTokenType::ALA},
      {-1,static_cast<resTokenType>(0)},
#line 80 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str105,    resTokenType::TPO},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 62 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str108,    resTokenType::DTH},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 54 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str113,    resTokenType::DHI},
      {-1,static_cast<resTokenType>(0)},
#line 39 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str115,    resTokenType::HIS},
      {-1,static_cast<resTokenType>(0)},
#line 98 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str117,     resTokenType::DU},
#line 119 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str118,    resTokenType::D3O},
      {-1,static_cast<resTokenType>(0)},
#line 37 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str120,    resTokenType::TYR},
      {-1,static_cast<resTokenType>(0)},
#line 44 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str122,    resTokenType::LYS},
#line 116 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str123,    resTokenType::H2O},
      {-1,static_cast<resTokenType>(0)},
#line 60 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str125,    resTokenType::DPR},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 105 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str128,    resTokenType::ACE},
      {-1,static_cast<resTokenType>(0)},
#line 123 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str130,    resTokenType::SPC},
#line 92 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str131,      resTokenType::N},
      {-1,static_cast<resTokenType>(0)},
#line 43 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str133,    resTokenType::GLN},
      {-1,static_cast<resTokenType>(0)},
#line 33 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str135,    resTokenType::CYS},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 38 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str138,    resTokenType::TRP},
      {-1,static_cast<resTokenType>(0)},
#line 51 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str140,    resTokenType::DCY},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 82 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str143,    resTokenType::MSE},
      {-1,static_cast<resTokenType>(0)},
#line 113 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str145,    resTokenType::SOL},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 45 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str148,    resTokenType::ARG},
      {-1,static_cast<resTokenType>(0)},
#line 50 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str150,    resTokenType::DAS},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 56 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str153,    resTokenType::DLE},
      {-1,static_cast<resTokenType>(0)},
#line 48 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str155,    resTokenType::DAR},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 40 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str158,    resTokenType::GLU},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 32 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str163,    resTokenType::THR},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 34 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str168,    resTokenType::MET},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 36 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str173,    resTokenType::PHE},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 99 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str177,     resTokenType::DN},
#line 58 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str178,    resTokenType::MED},
      {-1,static_cast<resTokenType>(0)},
#line 57 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str180,    resTokenType::DLY},
#line 90 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str181,      resTokenType::U},
      {-1,static_cast<resTokenType>(0)},
#line 108 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str183,    resTokenType::FMT},
      {-1,static_cast<resTokenType>(0)},
#line 26 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str185,    resTokenType::GLY},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 111 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str188,    resTokenType::NEH},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 30 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str193,    resTokenType::ILE},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 109 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str198,    resTokenType::E1H},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 66 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str203,    resTokenType::DNE},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 63 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str213,    resTokenType::DTR},
      {-1,static_cast<resTokenType>(0)},
#line 29 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str215,    resTokenType::LEU},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 106 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str218,    resTokenType::NH2},
      {-1,static_cast<resTokenType>(0)},
#line 64 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str220,    resTokenType::DTY},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 104 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str223,    resTokenType::NME},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 84 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str233,    resTokenType::PTR},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 35 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str238,    resTokenType::PRO},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 81 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str243,    resTokenType::SEC},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 79 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str253,    resTokenType::SEP},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 107 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str258,    resTokenType::FOR},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
      {-1,static_cast<resTokenType>(0)},
#line 31 "token.gperf"
      {(int)(size_t)&((struct stringpool_t *)0)->stringpool_str333,    resTokenType::SER}
    };

  if (len <= MAX_WORD_LENGTH && len >= MIN_WORD_LENGTH)
    {
      unsigned int key = hash (str, len);

      if (key <= MAX_HASH_VALUE)
        {
          int o = wordlist[key].string_offset;
          if (o >= 0)
            {
              const char *s = o + stringpool;

              if (*str == *s && !strncmp (str + 1, s + 1, len - 1) && s[len] == '\0')
                return &wordlist[key];
            }
        }
    }
  return 0;
}
#line 124 "token.gperf"

}

// #ifndef LOOKUP_IDENTIFIER_DEFINED
// #define LOOKUP_IDENTIFIER_DEFINED

inline resTokenType res_name_table(const char* res_name, std::size_t size) noexcept {
  const keywordEntry *entry = resNameKeyword::look_up(res_name, size);
  if (entry) {
    return entry->type;
  } else {
    return resTokenType::UNKNOWN;
  }
}
}
// #endif
#endif
