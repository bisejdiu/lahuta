#ifndef LAHUTA_DISTANCES_CONVERT_HPP
#define LAHUTA_DISTANCES_CONVERT_HPP

#include <stdexcept>
#include <type_traits>
#include <vector>

#include "rdkit/Geometry/point.h"

// clang-format off
namespace lahuta::dist {

inline std::vector<double> to_interleaved_xyz(const RDGeom::POINT3D_VECT &points) {
  std::vector<double> out;
  out.reserve(points.size() * 3);
  for (const auto &p : points) {
    out.push_back(p.x);
    out.push_back(p.y);
    out.push_back(p.z);
  }
  return out;
}

template <typename T>
std::vector<T> to_interleaved_xyz(const std::vector<std::vector<T>> &points) {
  static_assert(std::is_floating_point_v<T>, "T must be floating point");

  std::vector<T> out;
  out.reserve(points.size() * 3);
  for (const auto &row : points) {
    if (row.size() < 3) {
      throw std::invalid_argument("points must have at least 3 coordinates");
    }
    out.push_back(row[0]);
    out.push_back(row[1]);
    out.push_back(row[2]);
  }
  return out;
}

template std::vector<float>  to_interleaved_xyz<float> (const std::vector<std::vector<float>> &);
template std::vector<double> to_interleaved_xyz<double>(const std::vector<std::vector<double>> &);

}  // namespace lahuta::dist

#endif // LAHUTA_DISTANCES_CONVERT_HPP
