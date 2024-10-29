#ifndef LAHUTA_PISTACKING_HPP
#define LAHUTA_PISTACKING_HPP

#include "contacts/hydrogen_bonds.hpp"
#include "nn.hpp"

namespace lahuta {

class Luni;

void find_pistacking(const Luni *luni, GeometryOptions opts, Contacts &contacts);

} // namespace lahuta

#endif // LAHUTA_PISTACKING_HPP
