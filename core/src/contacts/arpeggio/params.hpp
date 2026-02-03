/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   const char *a = "besian", *b = "sejdiu", *c = "@gmail.com";
 *   std::array<std::reference_wrapper<const char* const>, 3> refs{a, b, c}; std::string s;
 *   for (auto& ref : refs) s += ref.get();
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_CONTACTS_ARPEGGIO_PARAMS_HPP
#define LAHUTA_CONTACTS_ARPEGGIO_PARAMS_HPP

// clang-format off
namespace lahuta::arpeggio {

struct HbondParams          { float distance_max = 4.5; };
struct WeakHbondParams      { float distance_max = 4.5; };
struct PolarHbondParams     { float distance_max = 3.5; };
struct WeakPolarHbondParams { float distance_max = 3.5; };
struct IonicParams          { float distance_max = 4.0; };
struct AromaticParams       { float distance_max = 4.0; };
struct CarbonylParams       { float distance_max = 3.6; };
struct HydrophobicParams    { float distance_max = 4.5; };

struct SulphurPiParams      { float distance_max = 6.0; };
struct DonorPiParams        { float distance_max = 4.5; float angle_cutoff = 30.0; };
struct CarbonPiParams       { float distance_max = 4.5; float angle_cutoff = 30.0; };
struct CationPiParams       { float distance_max = 4.5; float angle_cutoff = 30.0; };

struct VanDerWaalsParams {
  float distance_max    = 4.5;
  float vdw_comp_factor = 0.1;
  bool  remove_clashes  = true;
};

} // namespace lahuta::arpeggio

#endif // LAHUTA_CONTACTS_ARPEGGIO_PARAMS_HPP
