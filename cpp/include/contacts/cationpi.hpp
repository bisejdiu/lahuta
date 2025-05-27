#ifndef LAHUTA_CATIONPI_HPP
#define LAHUTA_CATIONPI_HPP

#include "entities/contact.hpp"

namespace lahuta {

class Topology;

struct CationPiParams {
  double distance_max = 6.0;
  double offset_max = 2.2;
};

ContactSet find_cationpi(const Topology& topology, const CationPiParams& params = CationPiParams{});

} // namespace lahuta

#endif // LAHUTA_CATIONPI_HPP
