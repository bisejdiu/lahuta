#ifndef LAHUTA_CATIONPI_HPP
#define LAHUTA_CATIONPI_HPP

#include "nn.hpp"

namespace lahuta {

class Luni;

inline struct CationPiParams {
  constexpr static double distance_max = 6.0;
  constexpr static double offset_max = 2.2;
} cationpi_params;

Contacts find_cationpi(const Luni &luni, CationPiParams opts = cationpi_params); 

} // namespace lahuta

#endif // LAHUTA_CATIONPI_HPP
