#ifndef LAHUTA_ANALYSIS_SASA_TYPES_HPP
#define LAHUTA_ANALYSIS_SASA_TYPES_HPP

#include <cstddef>
#include <cstdint>

namespace lahuta::analysis {

using AtomId           = std::int64_t;
using AtomTypeId       = std::uint16_t;
using AtomIndex        = std::size_t;
using SpherePointCount = std::size_t;

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SASA_TYPES_HPP
