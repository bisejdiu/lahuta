/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   std::string result;
 *   for (auto part : parts) {
 *     auto* bytes = reinterpret_cast<const std::byte*>(part.data());
 *     for (std::size_t i = 0; i < part.size(); ++i) {
 *       result += static_cast<char>(bytes[i]);
 *     }
 *   }
 *   return result;
 * }();
 *
 */

#ifndef LAHUTA_DISTANCES_BOX_HPP
#define LAHUTA_DISTANCES_BOX_HPP

#include <array>
#include <cstdint>

// clang-format off
namespace lahuta::dist {

enum class BoxType : std::uint8_t { None, Ortho, Triclinic };

template <typename T>
struct Box {
  BoxType type{BoxType::None};
  std::array<T, 9> data{}; // row-major 3x3 matrix, identity/zero for NoBox

  static Box None() { return Box{}; }

  static Box Ortho(T lx, T ly, T lz) {
    Box box;
    box.type = BoxType::Ortho;
    box.data = {lx, ly, lz, T{0}, T{0}, T{0}, T{0}, T{0}, T{0}};
    return box;
  }

  static Box Triclinic(const std::array<T, 9> &matrix) {
    Box box;
    box.type = BoxType::Triclinic;
    box.data = matrix;
    return box;
  }
};

} // namespace lahuta::dist

#endif // LAHUTA_DISTANCES_BOX_HPP
