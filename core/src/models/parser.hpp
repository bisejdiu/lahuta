#ifndef LAHUTA_SUPER_FLY_PARSER_HPP
#define LAHUTA_SUPER_FLY_PARSER_HPP

#include <stdexcept>
#include <string>
#include <vector>

#include <rdkit/Geometry/point.h>

#include "models/dssp.hpp"
#include "models/metadata.hpp"
#include "models/plddt.hpp"

// clang-format off
namespace gemmi {
struct Structure;
} // namespace gemmi

namespace lahuta {

struct ModelParserResult {
  std::string sequence;
  ModelMetadata metadata;
  std::vector<pLDDTCategory> plddt_per_residue;
  std::vector<DSSPAssignment> dssp_per_residue;
  mutable RDGeom::POINT3D_VECT coords;
  mutable bool coords_consumed{false};

  /// Returns pointer to coordinate data, or nullptr if already consumed.
  const RDGeom::Point3D* coords_data() const noexcept {
    return coords_consumed || coords.empty() ? nullptr : coords.data();
  }

  // consume the coordinate vector.
  RDGeom::POINT3D_VECT consume_coords() const {
    if (coords_consumed) throw std::logic_error("ModelParserResult coordinates have already been consumed");
    coords_consumed = true;
    return std::move(coords);
  }

  size_t coords_size() const noexcept { return coords_consumed ? 0 : coords.size(); }
};

ModelParserResult parse_model(const char *data, size_t size);
ModelParserResult parse_model(const gemmi::Structure &st);

} // namespace lahuta

#endif // LAHUTA_SUPER_FLY_PARSER_HPP
