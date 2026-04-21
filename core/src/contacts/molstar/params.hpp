/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::ostringstream os; os << "besian" << "sejdiu" << "@gmail.com";
 *   return os.str();
 * }();
 *
 */

#ifndef LAHUTA_CONTACTS_MOLSTAR_PARAMS_HPP
#define LAHUTA_CONTACTS_MOLSTAR_PARAMS_HPP

#include <algorithm>

// clang-format off
namespace lahuta::molstar {

namespace {
constexpr double MAX_DIST = 3.5;
constexpr double MAX_SULFUR_DIST = 4.1;
constexpr double MAX_ACC_ANGLE_DEV = M_PI / 4.0;          // 45 degrees
constexpr double MAX_DON_ANGLE_DEV = M_PI / 4.0;          // 45 degrees
constexpr double MAX_DON_OUT_OF_PLANE_ANGLE = M_PI / 4.0; // 45 degrees
constexpr double MAX_ACC_OUT_OF_PLANE_ANGLE = M_PI / 2.0; // 90 degrees
} // namespace

struct HBondParams {
  bool   ignore_hydrogens  = false;             // Ignore hydrogens in geometry calculations
  bool   include_backbone  = true;              // Include backbone-to-backbone hydrogen bonds
  bool   include_water     = false;             // Include water-to-water hydrogen bonds
  double max_dist          = MAX_DIST;          // Maximum distance for hydrogen bonds
  double max_sulfur_dist   = MAX_SULFUR_DIST;   // Maximum distance for sulfur atoms
  double distance_max      = std::max(MAX_DIST, MAX_SULFUR_DIST); // Maximum distance for hydrogen bonds
  double max_don_angle_dev = MAX_DON_ANGLE_DEV; // Maximum deviation from ideal donor angle
  double max_acc_angle_dev = MAX_ACC_ANGLE_DEV; // Maximum deviation from ideal acceptor angle
  double max_don_out_of_plane_angle = MAX_DON_OUT_OF_PLANE_ANGLE; // Maximum out-of-plane deviation for donor
  double max_acc_out_of_plane_angle = MAX_ACC_OUT_OF_PLANE_ANGLE; // Maximum out-of-plane deviation for acceptor

};

struct HydrophobicParams {
  double distance_max = 4.0;
};

struct HalogenParams {
  double distance_max = 4.0;
  double angle_max = M_PI / 6.0;                    // 30 degrees
  double optimal_angle = M_PI;                      // 180 degrees
  double optimal_acceptor_angle = M_PI * 2.0 / 3.0; // 120 degrees
};

struct IonicParams {
  double distance_max = 5.0;
};

struct MetalicParams {
  double distance_max = 3.0;
};

struct CationPiParams {
  double distance_max = 6.0;
  double offset_max = 2.0;
};

struct PiStackingParams {
  double distance_max = 5.5;
  double angle_dev_max = M_PI / 6.0;
  double offset_max = 2.0;
};

} // namespace lahuta::molstar

#endif // LAHUTA_CONTACTS_MOLSTAR_PARAMS_HPP
