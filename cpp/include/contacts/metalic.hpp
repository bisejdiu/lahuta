#ifndef LAHUTA_METALIC_HPP
#define LAHUTA_METALIC_HPP

#include "entities/contact.hpp"
#include <GraphMol/RWMol.h>

namespace lahuta {

class Topology;

struct MetalicParams {
  double distance_max = 3.0;
};

ContactSet find_metalic(const Topology &topology, const MetalicParams &params = MetalicParams{});

} // namespace lahuta

#endif // LAHUTA_METALIC_HPP
