#ifndef LAHUTA_SUPER_FLY_PARSER_HPP
#define LAHUTA_SUPER_FLY_PARSER_HPP

#include "Geometry/point.h"
#include <string>
#include <vector>

namespace lahuta {

// POINT3D_VECT is a vector of Point3D objects.
// class Point3D {
//  public:
//   double x{0.0};
//   double y{0.0};
//   double z{0.0};
// };

struct ModelParserResult {
  std::string sequence;
  RDGeom::POINT3D_VECT coords;

  std::vector<std::string> get_sequence() const {
    std::vector<std::string> result;
    result.reserve(sequence.size());
    for (char c : sequence) {
      result.push_back(std::string(1, c));
    }
    return result;
  }
};

ModelParserResult parse_model(const char *data, size_t size);

} // namespace lahuta

#endif // LAHUTA_SUPER_FLY_PARSER_HPP
