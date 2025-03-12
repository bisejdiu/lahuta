#ifndef LAHUTA_IONIC_HPP
#define LAHUTA_IONIC_HPP

#include "neighbors.hpp"

namespace lahuta {

class Luni;

inline struct IonicParams {
  constexpr static double distance_max = 5.0;
} ionic_params;

Contacts find_ionic(const Luni &luni, IonicParams opts = ionic_params);

} // namespace lahuta

#endif // LAHUTA_IONIC_HPP
