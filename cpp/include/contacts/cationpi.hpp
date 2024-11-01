#ifndef LAHUTA_CATIONPI_HPP
#define LAHUTA_CATIONPI_HPP

#include "nn.hpp"

namespace lahuta {

class Luni;

inline struct CationPiParams {
  constexpr static double distance_max = 6.0;
  constexpr static double offset_max = 2.2;
} cationpi_params;

RDGeom::Point3D project_on_plane(const RDGeom::Point3D &vector, const RDGeom::Point3D &plane_normal);

double compute_in_plane_offset(
    const RDGeom::Point3D &pos_a, const RDGeom::Point3D &pos_b, const RDGeom::Point3D &normal);

Contacts find_cationpi(const Luni &luni, CationPiParams opts = cationpi_params); 

} // namespace lahuta

#endif // LAHUTA_CATIONPI_HPP
