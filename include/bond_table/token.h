#ifndef TOKENS_H
#define TOKENS_H

#include <cstdint>

namespace lahuta {

enum class amino_acid_names_l : std::uint8_t {
  HIS,  // histidine
  ARG,  // arginine
  LYS,  // lysine
  ILE,  // isoleucine
  PHE,  // phenylalanine
  LEU,  // leucine
  TRP,  // tryptophan
  ALA,  // alanine
  MET,  // methionine
  PRO,  // proline
  CYS,  // cysteine
  ASN,  // asparagine
  VAL,  // valine
  GLY,  // glycine
  SER,  // serine
  GLN,  // glutamine
  TYR,  // tyrosine
  ASP,  // aspartic acid
  GLU,  // glutamic acid
  THR,  // threonine
  SEC,  // selenocysteine
  PYL,  // pyrrolysine
  UNK,  // unknown amino acid from CCD
  MSE,  // selenomethionine
  SEP,  // phosphoserine
  TPO,  // phosphothreonine
  PTR,  // phosphotyrosine
  PCA,  // pyrrolidone carboxylic acid
  HYP,  // hydroxyproline
  HSD,  // histidine delta tautomer
  HSE,  // histidine epsilon tautomer
  HSP,  // histidine protonated
  LSN,  // lysine protonated
  ASPP, // aspartic acid protonated
  GLUP, // glutamic acid protonated
  HID,  // histidine delta tautomer
  HIE,  // histidine epsilon tautomer
  HIP,  // histidine protonated on both nitrogens
  LYN,  // lysine protonated
  ASH,  // aspartic acid protonated
  GLH,  // glutamic acid protonated
};

enum class amino_acid_names_d : std::uint8_t {
  DAL, // D-Alanine
  DAR, // D-Arginine
  DSG, // D-Asparagine
  DAS, // D-Aspartic Acid
  DCY, // D-Cysteine
  DGL, // D-Glutamic Acid
  DGN, // D-Glutamine
  DHI, // D-Histidine
  DIL, // D-Isoleucine
  DLE, // D-Leucine
  DLY, // D-Lysine
  MED, // D-Methionine
  DPN, // D-Phenylalanine
  DPR, // D-Proline
  DSN, // D-Serine
  DTH, // D-Threonine
  DTR, // D-Tryptophan
  DTY, // D-Tyrosine
  DVA, // D-Valine
  DNE, // D-Norleucine
};

enum class common_protein_caps : std::uint8_t {
  NME, // N-terminal methyl
  ACE, // N-terminal acetyl
  NH2, // C-terminal amide
  FOR, // formyl
  FMT, // formate
  // E1H, // GFP backbone fragmentation in 2G16
  // HOA, // complexes zinc
  // NEH, // ubiquitine linker
  // MOH, // part of peptidomimetics
};

enum class resTokenType : std::uint8_t {
  // Amino Acid Names L
  GLY,  // glycine
  ALA,  // alanine
  VAL,  // valine
  LEU,  // leucine
  ILE,  // isoleucine
  SER,  // serine
  THR,  // threonine
  CYS,  // cysteine
  MET,  // methionine
  PRO,  // proline
  PHE,  // phenylalanine
  TYR,  // tyrosine
  TRP,  // tryptophan
  HIS,  // histidine
  GLU,  // glutamic acid
  ASP,  // aspartic acid
  ASN,  // asparagine
  GLN,  // glutamine
  LYS,  // lysine
  ARG,  // arginine
  PCA,  // pyrrolidone carboxylic acid

  // Amino Acid Names D
  DAL, // D-Alanine
  DAR, // D-Arginine
  DSG, // D-Asparagine
  DAS, // D-Aspartic Acid
  DCY, // D-Cysteine
  DGL, // D-Glutamic Acid
  DGN, // D-Glutamine
  DHI, // D-Histidine
  DIL, // D-Isoleucine
  DLE, // D-Leucine
  DLY, // D-Lysine
  MED, // D-Methionine
  DPN, // D-Phenylalanine
  DPR, // D-Proline
  DSN, // D-Serine
  DTH, // D-Threonine
  DTR, // D-Tryptophan
  DTY, // D-Tyrosine
  DVA, // D-Valine
  DNE, // D-Norleucine

  // Charmm-FF and Amber-FF Residue Names
  HSD,  // histidine delta tautomer
  HSE,  // histidine epsilon tautomer
  HSP,  // histidine protonated
  HID,  // histidine delta tautomer
  HIE,  // histidine epsilon tautomer
  HIP,  // histidine protonated on both nitrogens
  GLUP, // glutamic acid protonated
  ASH,  // aspartic acid protonated
  ASPP, // aspartic acid protonated
  GLH,  // glutamic acid protonated
  LSN,  // lysine protonated
  LYN,  // lysine protonated

  // Modified Residues
  SEP,  // phosphoserine
  TPO,  // phosphothreonine
  SEC,  // selenocysteine
  MSE,  // selenomethionine
  HYP,  // hydroxyproline
  PTR,  // phosphotyrosine
  PYL,  // pyrrolysine

  // RNA Base Names
  A,    // adenine
  C,    // cytosine
  G,    // guanine
  T,    // thymine
  U,    // uracil
  I,    // inosine
  N,    // nucleotide ?

  // DNA Base Names
  DA,   // deoxyadenosine
  DC,   // deoxycytidine
  DG,   // deoxyguanosine
  DT,   // deoxythymidine
  DI,   // deoxyinosine
  DU,   // deoxyuridine
  DN,   // deoxynucleotide

  // Peptide Base Names
  APN,  // adenosine phosphate
  CPN,  // cytidine phosphate
  TPN,  // thymidine phosphate
  GPN,  // guanosine phosphate

  // Common Protein Caps
  NME,  // N-terminal methyl
  ACE,  // N-terminal acetyl
  NH2,  // C-terminal amide
  FOR,  // formyl
  FMT,  // formate
  E1H,  // GFP backbone fragmentation in 2G16
  HOA,  // complexes zinc
  NEH,  // ubiquitine linker
  MOH,  // part of peptidomimetics

  // Water Names
  SOL,
  WAT,
  HOH,
  H2O,
  W,
  DOD,
  D3O,
  TIP,
  TIP3,
  TIP4,
  SPC,

  UNKNOWN,
};


// #define RES_KEYWORD_TOKENS_TABLE 	\
// 	{"GLY", resTokenType::GLY},	\
// 	{"ALA", resTokenType::ALA},	\
// 	{"VAL", resTokenType::VAL},	\
// 	{"LEU", resTokenType::LEU},	\
// 	{"ILE", resTokenType::ILE},	\
// 	{"SER", resTokenType::SER},	\
// 	{"THR", resTokenType::THR},	\
// 	{"CYS", resTokenType::CYS},	\
// 	{"MET", resTokenType::MET},	\
// 	{"PRO", resTokenType::PRO},	\
// 	{"PHE", resTokenType::PHE},	\
// 	{"TYR", resTokenType::TYR},	\
// 	{"TRP", resTokenType::TRP},	\
// 	{"HIS", resTokenType::HIS},	\
// 	{"GLU", resTokenType::GLU},	\
// 	{"ASP", resTokenType::ASP},	\
// 	{"ASN", resTokenType::ASN},	\
// 	{"GLN", resTokenType::GLN},	\
// 	{"LYS", resTokenType::LYS},	\
// 	{"ARG", resTokenType::ARG},
//
// struct resKeywordToken {
//   std::string_view keyword;
//   resTokenType token;
// };
//
// constexpr resKeywordToken resTokensTable[] = {
//   RES_KEYWORD_TOKENS_TABLE
// };

typedef resTokenType res_name_table_fn(const char *res_name,
                                         std::size_t size) noexcept;

extern resTokenType res_name_table(const char *res_name,
                                     std::size_t size) noexcept;

} // namespace lahuata

#endif
