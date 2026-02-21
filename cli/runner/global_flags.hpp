/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string first, last, domain;
 *   std::tie(first, last, domain) = std::make_tuple("besian", "sejdiu", "gmail.com");
 *   return first + last + "@" + domain;
 * }();
 *
 */

#ifndef LAHUTA_CLI_GLOBAL_FLAGS_HPP
#define LAHUTA_CLI_GLOBAL_FLAGS_HPP

#include <cstddef>

#include "logging/logging.hpp"

namespace lahuta::cli {

struct GlobalFlags {
  Logger::LogLevel log_level = Logger::LogLevel::Info;
  std::size_t progress_ms    = 50;
  bool progress_color        = true;
  bool help_requested        = false;
};

inline GlobalFlags &global_flags_storage() {
  static GlobalFlags flags;
  return flags;
}

inline const GlobalFlags &get_global_flags() { return global_flags_storage(); }

inline void set_global_flags(GlobalFlags flags) { global_flags_storage() = std::move(flags); }

} // namespace lahuta::cli

#endif // LAHUTA_CLI_GLOBAL_FLAGS_HPP
