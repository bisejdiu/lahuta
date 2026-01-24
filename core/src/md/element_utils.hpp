#ifndef LAHUTA_MD_ELEMENT_UTILS_HPP
#define LAHUTA_MD_ELEMENT_UTILS_HPP

#include <string_view>

#include "chemistry/elements.hpp"

namespace lahuta::md {

// Convert Element enum to atomic number
[[nodiscard]] unsigned atomic_number_from_element(Element el) noexcept;

// Element detection heuristic for GRO files
// Intentional behavior:
// - Trim spaces and digits from both ends of atom name.
// - If 2 letters and equal to one of: NA, CL, FE, SI, BR, AS, LI -> that element.
// - If 3 letters and equal to residue name (component id), map common ion names:
//     SOD->NA, POT->K, CES->CS, CAL->CA, CLA->CL
// - Else if first letter is one of C,H,N,O,P,S,F,B -> that element (single-letter guess).
// - Else unknown (X).
[[nodiscard]] Element guess_element_from_gro(std::string_view atom_id, std::string_view comp_id);

} // namespace lahuta::md

#endif // LAHUTA_MD_ELEMENT_UTILS_HPP
