#ifndef LAHUTA_METALIC_HPP
#define LAHUTA_METALIC_HPP

#include "neighbors.hpp"
namespace lahuta {

struct MetalicParams {
  double distance_max = 3.0;
};

Contacts find_metalic(const Luni &luni, std::optional<MetalicParams> params = std::nullopt);

} // namespace lahuta

#endif // LAHUTA_METALIC_HPP
