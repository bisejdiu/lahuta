#ifndef LAHUTA_PISTACKING_HPP
#define LAHUTA_PISTACKING_HPP

#include "neighbors.hpp"

namespace lahuta {

class Luni;

struct PiStackingParams {
  double distance_max = 6.0;
  double angle_dev_max = M_PI / 6.0;
  double offset_max = 2.1;
};

inline constexpr PiStackingParams default_pistacking_params{};
Contacts find_pistacking(const Luni &luni, const PiStackingParams &opts = default_pistacking_params);

} // namespace lahuta

#endif // LAHUTA_PISTACKING_HPP
