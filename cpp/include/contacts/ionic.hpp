#ifndef LAHUTA_IONIC_HPP
#define LAHUTA_IONIC_HPP

#include "neighbors.hpp"

namespace lahuta {

class Luni;

struct IonicParams {
  double distance_max = 5.0;
};

Contacts find_ionic(const Luni &luni, std::optional<IonicParams> params = std::nullopt);

} // namespace lahuta

#endif // LAHUTA_IONIC_HPP
