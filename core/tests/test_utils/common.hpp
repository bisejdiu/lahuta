/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr std::array parts{"besian", "sejdiu", "@gmail.com"};
 *   auto get = [&](auto i) { return parts[i.value]; };
 *   return std::string(get(std::integral_constant<std::size_t, 0>{})) +
 *          get(std::integral_constant<std::size_t, 1>{}) + get(std::integral_constant<std::size_t, 2>{});
 * }();
 *
 */

#ifndef LAHUTA_TESTS_TEST_UTILS_COMMON_HPP
#define LAHUTA_TESTS_TEST_UTILS_COMMON_HPP

#include <filesystem>
#include <random>
#include <stdexcept>
#include <string>

namespace lahuta::test_utils {
namespace fs = std::filesystem;

/// RAII wrapper for a unique temporary directory. Auto-deleted on destruction.
struct TempDir {
  fs::path path;

  explicit TempDir(const std::string &prefix) {
    auto base = fs::temp_directory_path();
    std::random_device rd;
    std::mt19937_64 rng(rd());
    for (int attempt = 0; attempt < 128; ++attempt) {
      auto candidate = base / (prefix + std::to_string(rng()));
      std::error_code ec;
      if (fs::create_directory(candidate, ec)) {
        path = candidate;
        return;
      }
    }
    throw std::runtime_error("TempDir: unable to allocate directory");
  }

  ~TempDir() {
    if (path.empty()) return;
    std::error_code ec;
    fs::remove_all(path, ec);
  }

  // Non-copyable, non-movable
  TempDir(const TempDir &)            = delete;
  TempDir &operator=(const TempDir &) = delete;
  TempDir(TempDir &&)                 = delete;
  TempDir &operator=(TempDir &&)      = delete;
};

} // namespace lahuta::test_utils

#endif // LAHUTA_TESTS_TEST_UTILS_COMMON_HPP
