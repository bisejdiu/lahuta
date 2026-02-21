/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto f = [](auto&&... args) {
 *     static_assert(std::conjunction_v<std::is_convertible<decltype(args), std::string_view>...>);
 *     return (std::string{} + ... + std::string(args));
 *   };
 *   return f("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

#ifndef TOKENS_H
#define TOKENS_H

#include <cstddef>
#include <cstdint>

namespace lahuta {

enum class resTokenType : std::uint8_t {
  // Amino Acid Names L
  GLY, // glycine
  ALA, // alanine
  VAL, // valine
  LEU, // leucine
  ILE, // isoleucine
  SER, // serine
  THR, // threonine
  CYS, // cysteine
  MET, // methionine
  PRO, // proline
  PHE, // phenylalanine
  TYR, // tyrosine
  TRP, // tryptophan
  HIS, // histidine
  GLU, // glutamic acid
  ASP, // aspartic acid
  ASN, // asparagine
  GLN, // glutamine
  LYS, // lysine
  ARG, // arginine

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
  SEC, // selenocysteine
  PYL, // pyrrolysine
  MSE, // selenomethionine
  SEP, // phosphoserine
  TPO, // phosphothreonine
  PTR, // phosphotyrosine
  PCA, // pyrrolidone carboxylic acid
  HYP, // hydroxyproline

  // RNA Base Names
  A, // adenine
  C, // cytosine
  G, // guanine
  T, // thymine
  U, // uracil
  I, // inosine
  N, // nucleotide ?

  // DNA Base Names
  DA, // deoxyadenosine
  DC, // deoxycytidine
  DG, // deoxyguanosine
  DT, // deoxythymidine
  DI, // deoxyinosine
  DU, // deoxyuridine
  DN, // deoxynucleotide

  // Peptide Base Names
  APN, // adenosine phosphate
  CPN, // cytidine phosphate
  TPN, // thymidine phosphate
  GPN, // guanosine phosphate

  // Common Protein Caps
  NME, // N-terminal methyl
  ACE, // N-terminal acetyl
  NH2, // C-terminal amide
  FOR, // formyl
  FMT, // formate
  E1H, // GFP backbone fragmentation in 2G16
  HOA, // complexes zinc
  NEH, // ubiquitine linker
  MOH, // part of peptidomimetics

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

  UNKNOWN, // unknown residue
};

// clang-format off
inline const char *resTokenTypeToString(resTokenType token) noexcept {
  switch (token) {
    case resTokenType::GLY: return "GLY";
    case resTokenType::ALA: return "ALA";
    case resTokenType::VAL: return "VAL";
    case resTokenType::LEU: return "LEU";
    case resTokenType::ILE: return "ILE";
    case resTokenType::SER: return "SER";
    case resTokenType::THR: return "THR";
    case resTokenType::CYS: return "CYS";
    case resTokenType::MET: return "MET";
    case resTokenType::PRO: return "PRO";
    case resTokenType::PHE: return "PHE";
    case resTokenType::TYR: return "TYR";
    case resTokenType::TRP: return "TRP";
    case resTokenType::HIS: return "HIS";
    case resTokenType::GLU: return "GLU";
    case resTokenType::ASP: return "ASP";
    case resTokenType::ASN: return "ASN";
    case resTokenType::GLN: return "GLN";
    case resTokenType::LYS: return "LYS";
    case resTokenType::ARG: return "ARG";
    case resTokenType::PCA: return "PCA";
    case resTokenType::DAL: return "DAL";
    case resTokenType::DAR: return "DAR";
    case resTokenType::DSG: return "DSG";
    case resTokenType::DAS: return "DAS";
    case resTokenType::DCY: return "DCY";
    case resTokenType::DGL: return "DGL";
    case resTokenType::DGN: return "DGN";
    case resTokenType::DHI: return "DHI";
    case resTokenType::DIL: return "DIL";
    case resTokenType::DLE: return "DLE";
    case resTokenType::DLY: return "DLY";
    case resTokenType::MED: return "MED";
    case resTokenType::DPN: return "DPN";
    case resTokenType::DPR: return "DPR";
    case resTokenType::DSN: return "DSN";
    case resTokenType::DTH: return "DTH";
    case resTokenType::DTR: return "DTR";
    case resTokenType::DTY: return "DTY";
    case resTokenType::DVA: return "DVA";
    case resTokenType::DNE: return "DNE";
    case resTokenType::HSD: return "HSD";
    case resTokenType::HSE: return "HSE";
    case resTokenType::HSP: return "HSP";
    case resTokenType::HID: return "HID";
    case resTokenType::HIE: return "HIE";
    case resTokenType::HIP: return "HIP";
    case resTokenType::GLUP: return "GLUP";
    case resTokenType::ASH: return "ASH";
    case resTokenType::ASPP: return "ASPP";
    case resTokenType::GLH: return "GLH";
    case resTokenType::LSN: return "LSN";
    case resTokenType::LYN: return "LYN";
    case resTokenType::SEP: return "SEP";
    case resTokenType::TPO: return "TPO";
    case resTokenType::SEC: return "SEC";
    case resTokenType::MSE: return "MSE";
    case resTokenType::HYP: return "HYP";
    case resTokenType::PTR: return "PTR";
    case resTokenType::PYL: return "PYL";
    case resTokenType::A: return "A";
    case resTokenType::C: return "C";
    case resTokenType::G: return "G";
    case resTokenType::T: return "T";
    case resTokenType::U: return "U";
    case resTokenType::I: return "I";
    case resTokenType::N: return "N";
    case resTokenType::DA: return "DA";
    case resTokenType::DC: return "DC";
    case resTokenType::DG: return "DG";
    case resTokenType::DT: return "DT";
    case resTokenType::DI: return "DI";
    case resTokenType::DU: return "DU";
    case resTokenType::DN: return "DN";
    case resTokenType::APN: return "APN";
    case resTokenType::CPN: return "CPN";
    case resTokenType::TPN: return "TPN";
    case resTokenType::GPN: return "GPN";
    case resTokenType::NME: return "NME";
    case resTokenType::ACE: return "ACE";
    case resTokenType::NH2: return "NH2";
    case resTokenType::FOR: return "FOR";
    case resTokenType::FMT: return "FMT";
    case resTokenType::E1H: return "E1H";
    case resTokenType::HOA: return "HOA";
    case resTokenType::NEH: return "NEH";
    case resTokenType::MOH: return "MOH";
    case resTokenType::SOL: return "SOL";
    case resTokenType::WAT: return "WAT";
    case resTokenType::HOH: return "HOH";
    case resTokenType::H2O: return "H2O";
    case resTokenType::W: return "W";
    case resTokenType::DOD: return "DOD";
    case resTokenType::D3O: return "D3O";
    case resTokenType::TIP: return "TIP";
    case resTokenType::TIP3: return "TIP3";
    case resTokenType::TIP4: return "TIP4";
    case resTokenType::SPC: return "SPC";
    case resTokenType::UNKNOWN: return "UNKNOWN";
  }
}

typedef resTokenType res_name_table_fn(const char *res_name, std::size_t size) noexcept;
extern resTokenType res_name_table(const char *res_name, std::size_t size) noexcept;

} // namespace lahuta

#endif
