#ifndef LAHUTA_CATIONPI_HPP
#define LAHUTA_CATIONPI_HPP

#include "neighbors.hpp"

namespace lahuta {

class Luni;

struct CationPiParams {
  double distance_max = 6.0;
  double offset_max = 2.2;
};

Contacts find_cationpi(const Luni &luni, std::optional<CationPiParams> params = std::nullopt);

} // namespace lahuta

#endif // LAHUTA_CATIONPI_HPP
