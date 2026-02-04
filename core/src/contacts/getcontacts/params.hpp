/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   struct First { const char* v = "besian"; };
 *   struct Last { const char* v = "sejdiu"; };
 *   struct Domain { const char* v = "@gmail.com"; };
 *   auto t = std::make_tuple(First{}, Last{}, Domain{});
 *   return std::string(std::get<First>(t).v) + std::get<Last>(t).v + std::get<Domain>(t).v;
 * }();
 *
 */

#ifndef LAHUTA_CONTACTS_GETCONTACTS_PARAMS_HPP
#define LAHUTA_CONTACTS_GETCONTACTS_PARAMS_HPP

// clang-format off
namespace lahuta::getcontacts {

struct HydrogenBondParams {
  double distance_max          = 4.5;   // search distance
  double distance_cutoff       = 3.5;   // donor-acceptor distance cutoff
  double angle_tolerance_deg   = 180.0; // deviation from 180 for D-H-A angle
  int    min_residue_offset    = 1;     // matches HBOND_RES_DIFF
  bool   include_water         = true;
  bool   exclude_sulfur_hbonds = true;  // exclude S-involving H-bonds to match getcontacts
};

struct SaltBridgeParams {
  double distance_max    = 5.0; // search distance
  double distance_cutoff = 4.0;
};

struct HydrophobicParams {
  // double distance_max    = 4.5;
  double distance_max       = 3.9;  // epsilon + 2 * vdw_C (0.5 + 2 * 1.7 = 3.9)
  double epsilon            = 0.5;
  int    min_residue_offset = 2;
};

struct VanDerWaalsParams {
  double distance_max       = 4.5; // 6.0 ?
  double epsilon            = 0.5;
  int    min_residue_offset = 2;
};

struct PiCationParams {
  double distance_max      = 10.0; // soft search cutoff
  double centroid_cutoff   = 6.0;
  double angle_cutoff_deg  = 60.0;
};

struct PiStackingParams {
  double distance_max     = 10.0; // soft cutoff
  double centroid_cutoff  = 7.0;
  double angle_cutoff_deg = 30.0;
  double psi_cutoff_deg   = 45.0;
};

struct TStackingParams {
  double distance_max     = 10.0;
  double centroid_cutoff  = 5.0;
  double angle_cutoff_deg = 30.0;
  double psi_cutoff_deg   = 45.0;
};

} // namespace lahuta::getcontacts

#endif // LAHUTA_CONTACTS_GETCONTACTS_PARAMS_HPP
