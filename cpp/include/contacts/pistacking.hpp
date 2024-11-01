#ifndef LAHUTA_PISTACKING_HPP
#define LAHUTA_PISTACKING_HPP

#include "contacts/hydrogen_bonds.hpp"
#include "nn.hpp"

namespace lahuta {

class Luni;

constexpr double pistacking_dist_max = 5.5; // Maximum distance for π-stacking interactions
constexpr double offset_max = 2.1;          // Maximum offset for π-stacking interactions
constexpr double AngleDevMax = M_PI / 6.0;  // 30 degrees in radians

void find_pistacking(const Luni *luni, GeometryOptions opts, Contacts &contacts);

} // namespace lahuta

#endif // LAHUTA_PISTACKING_HPP
