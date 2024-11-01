#ifndef LAHUTA_METALIC_HPP
#define LAHUTA_METALIC_HPP

#include "nn.hpp"
namespace lahuta {

bool is_metalic(AtomType at1, AtomType at2);

inline struct MetalicParams {
  constexpr static double distance_max = 3.0;
} metalic_params;

Contacts find_metalic(const Luni &luni, MetalicParams opts = metalic_params);

} // namespace lahuta

#endif // LAHUTA_METALIC_HPP
