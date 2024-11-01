#ifndef LAHUTA_METALIC_HPP
#define LAHUTA_METALIC_HPP

#include "nn.hpp"
#include "contacts/hydrogen_bonds.hpp"
namespace lahuta {

bool is_metalic(AtomType at1, AtomType at2);

void find_metalic(const Luni *luni, GeometryOptions opts, Contacts &contacts);

} // namespace lahuta

#endif // LAHUTA_METALIC_HPP
