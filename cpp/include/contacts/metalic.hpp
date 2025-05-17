#ifndef LAHUTA_METALIC_HPP
#define LAHUTA_METALIC_HPP

#include "neighbors.hpp"
namespace lahuta {

inline struct MetalicParams {
  constexpr static double distance_max = 3.0;
} metalic_params;

Contacts find_metalic(const Luni &luni, MetalicParams opts = metalic_params);

} // namespace lahuta

#endif // LAHUTA_METALIC_HPP
