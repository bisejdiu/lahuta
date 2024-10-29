#ifndef LAHUTA_IONIC_HPP
#define LAHUTA_IONIC_HPP

#include "nn.hpp"
#include "contacts/hydrogen_bonds.hpp"

namespace lahuta {

class Luni;

void find_ionic(const Luni &luni, GeometryOptions opts, Contacts &container);

} // namespace lahuta

#endif // LAHUTA_IONIC_HPP
