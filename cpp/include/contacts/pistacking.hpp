#ifndef LAHUTA_PISTACKING_HPP
#define LAHUTA_PISTACKING_HPP

#include "nn.hpp"

namespace lahuta {

class Luni;

inline struct PiStackingParams {
  constexpr static double distance_max = 6.0;
  constexpr static double angle_dev_max = M_PI / 6.0;
  constexpr static double offset_max = 2.1;
} pistacking_params;

Contacts find_pistacking(const Luni &luni, PiStackingParams opts = pistacking_params);

} // namespace lahuta

#endif // LAHUTA_PISTACKING_HPP
