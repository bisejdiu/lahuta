#ifndef LAHUTA_IONIC_HPP
#define LAHUTA_IONIC_HPP

#include "entities/contact.hpp"

namespace lahuta {

class Luni;
class Topology;

struct IonicParams {
  double distance_max = 5.0;
};

ContactSet find_ionic(const Topology& topology, const IonicParams& params = IonicParams{});

} // namespace lahuta

#endif // LAHUTA_IONIC_HPP
