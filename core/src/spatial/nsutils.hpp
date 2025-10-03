#ifndef LAHUTA_NSUTILS_HPP
#define LAHUTA_NSUTILS_HPP

#include <array>
#include <cstddef>

#include <rdkit/GraphMol/MonomerInfo.h>

namespace lahuta::ns_utils {

// Hard cap to avoid excessive pre-allocation
constexpr std::size_t PAIRS_RESERVE_CAP = 5'000'000;

// Estimate the number of undirected neighbor pairs to pre-reserve for self-search (i<j)
// Uses uniform density assumption. Intentionally conservative, with hard cap.
std::size_t estimate_pairs(std::size_t N, float cutoff, const std::array<float, 3> &box);

// Estimate the number of directed neighbor pairs to pre-reserve for cross-search (m queries against N points)
std::size_t estimate_pairs(std::size_t M, std::size_t N, float cutoff, const std::array<float, 3> &box);

} // namespace lahuta::ns_utils

#endif // LAHUTA_NSUTILS_HPP
