#ifndef LAHUTA_VERSION_HPP
#define LAHUTA_VERSION_HPP

#include <string_view>

#ifndef LAHUTA_VERSION_STRING
#define LAHUTA_VERSION_STRING "0.0.0"
#endif

#ifndef LAHUTA_VERSION_SUFFIX
#define LAHUTA_VERSION_SUFFIX ""
#endif

#ifndef LAHUTA_VERSION_MAJOR
#define LAHUTA_VERSION_MAJOR 0
#endif

#ifndef LAHUTA_VERSION_MINOR
#define LAHUTA_VERSION_MINOR 0
#endif

#ifndef LAHUTA_VERSION_PATCH
#define LAHUTA_VERSION_PATCH 0
#endif

namespace lahuta {
inline constexpr std::string_view version{LAHUTA_VERSION_STRING};
inline constexpr int version_major{LAHUTA_VERSION_MAJOR};
inline constexpr int version_minor{LAHUTA_VERSION_MINOR};
inline constexpr int version_patch{LAHUTA_VERSION_PATCH};
inline constexpr std::string_view version_suffix{LAHUTA_VERSION_SUFFIX};
inline constexpr bool version_is_prerelease{version_suffix.size() > 0};
} // namespace lahuta

#endif // LAHUTA_VERSION_HPP
