#ifndef LAHUTA_CATIONPI_HPP
#define LAHUTA_CATIONPI_HPP

#include "contacts/hydrogen_bonds.hpp"
#include "nn.hpp"
namespace lahuta {

class Luni;

RDGeom::Point3D project_on_plane(const RDGeom::Point3D &vector, const RDGeom::Point3D &plane_normal);

double compute_in_plane_offset(
    const RDGeom::Point3D &pos_a, const RDGeom::Point3D &pos_b, const RDGeom::Point3D &normal);

void find_cationpi(const Luni *luni, GeometryOptions opts, Contacts &contacts);

} // namespace lahuta

#endif // LAHUTA_CATIONPI_HPP
