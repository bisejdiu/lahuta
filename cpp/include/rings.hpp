#ifndef LAHUTA_RINGS_HPP
#define LAHUTA_RINGS_HPP

#include "GraphMol/MonomerInfo.h"
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/RWMol.h>
#include <vector>

namespace lahuta {

// TODO: 1. Store the ring type (aromatic, aliphatic, etc.)

enum class RingGroup { NONE, AROMATIC, ALIPHATIC };
enum class RingType { None };

struct RingData {
  /*std::vector<int> atom_ids;*/
  std::vector<const RDKit::Atom *> atoms;
  RDGeom::Point3D center;
  RDGeom::Point3D norm;

  size_t get_id() const { return id; }

private:
  size_t id;

public:
  RingData() = default;
  /*RingData(RDGeom::Point3D center_, RDGeom::Point3D norm_, std::vector<int>
   * atom_ids_)*/
  /*    : center(center_), norm(norm_), atom_ids(atom_ids_) {}*/

  explicit RingData(RDGeom::Point3D center_, RDGeom::Point3D norm_, std::vector<const RDKit::Atom *> atoms_)
      : center(center_), norm(norm_), atoms(atoms_) {}

  explicit RingData(RDGeom::Point3D center_, RDGeom::Point3D norm_, std::vector<const RDKit::Atom *> atoms_, size_t id_)
      : center(center_), norm(norm_), atoms(atoms_), id(id_) {}

  std::vector<int> atom_ids() const {
    std::vector<int> ids;
    ids.reserve(atoms.size());
    for (const auto *atom : atoms) {
      ids.push_back(atom->getIdx());
    }
    return ids;
  }

  /// Compute the angle between the provided point and the ring
  double compute_angle(const std::vector<double> &_point) const {
    RDGeom::Point3D point(_point[0], _point[1], _point[2]);
    auto vector_point_to_plane = point - center;
    vector_point_to_plane.normalize();

    double cos_theta = vector_point_to_plane.dotProduct(norm);
    cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

    double theta_radians = std::acos(cos_theta); // in radians
    return theta_radians * (180.0 / M_PI);
  }

  // Function to compute the angle between the point and the plane
  double compute_angle_(const std::vector<double> _point) const {
    RDGeom::Point3D point(_point[0], _point[1], _point[2]);
    RDGeom::Point3D vector_point_to_plane = point - center;

    RDGeom::Point3D normalized_line_direction = vector_point_to_plane;
    normalized_line_direction.normalize();
    double dot_product = normalized_line_direction.dotProduct(norm);

    double raw_angle = std::acos(dot_product);
    double adjusted_angle = std::copysign(raw_angle, dot_product);
    if (adjusted_angle < 0.0) {
      adjusted_angle += M_PI;
    }
    return adjusted_angle * (180.0 / M_PI);
  }
};

struct RingDataVec {
  std::vector<RingData> rings;

  const std::vector<RingData> &get_data() const { return rings; }

  RingData &operator[](size_t index) { return rings[index]; }
  const RingData &operator[](size_t index) const { return rings[index]; }
  int size() const { return rings.size(); }

  RingDataVec() = default;
  RingDataVec(RingDataVec &&other) noexcept : rings(std::move(other.rings)) {}
  RingDataVec &operator=(RingDataVec &&other) noexcept {
    if (this != &other) {
      rings = std::move(other.rings);
    }
    return *this;
  }

  RingDataVec(const RingDataVec &other) : rings(other.rings) {}
  RingDataVec &operator=(const RingDataVec &other) {
    if (this != &other) {
      rings = other.rings;
    }
    return *this;
  }

  const RDGeom::POINT3D_VECT centers() const {
    RDGeom::POINT3D_VECT centers;
    for (const auto &ring : rings) {
      centers.push_back(ring.center);
    }
    return centers;
  }

  const RDGeom::POINT3D_VECT positions() const { return centers(); }

  std::vector<std::vector<double>> norm1() const {
    std::vector<std::vector<double>> normals;
    for (const auto &ring : rings) {
      normals.push_back({ring.norm.x, ring.norm.y, ring.norm.z});
    }
    return normals;
  }

  /*std::vector<int> root_atom_ids() const {*/
  /*  std::vector<int> root_atom_ids;*/
  /*  for (const auto &ring : rings) {*/
  /*    root_atom_ids.push_back(ring.get_root_atom_id());*/
  /*  }*/
  /*  return root_atom_ids;*/
  /*}*/

  std::vector<double>
  compute_angles(const std::vector<int> &ring_indices, const std::vector<std::vector<double>> &points) const {
    std::vector<double> angles;
    angles.reserve(ring_indices.size());
    for (size_t i = 0; i < ring_indices.size(); ++i) {
      angles.push_back(rings[ring_indices[i]].compute_angle(points[i]));
    }
    return angles;
  }
};

class RingProps {
public:
  static void
  compute_center(const RDKit::RWMol *mol, const std::vector<int> &ring_atom_ids, RDGeom::Point3D &center) {
    const RDKit::Conformer &conf = mol->getConformer();
    center = RDGeom::Point3D{0.0, 0.0, 0.0};
    for (int idx : ring_atom_ids) {
      center += conf.getAtomPos(idx);
    }
    center /= static_cast<double>(ring_atom_ids.size());
  }

  static void compute_normal(
      const RDKit::RWMol *mol, const std::vector<int> &ring_atom_ids, const RDGeom::Point3D &center,
      RDGeom::Point3D &norm) {
    if (ring_atom_ids.size() < 3) {
      throw std::invalid_argument("Ring must contain at least 3 atoms to compute a normal.");
    }

    const RDKit::Conformer &conf = mol->getConformer();

    for (size_t i = 0; i < ring_atom_ids.size(); ++i) {
      RDGeom::Point3D v1 = conf.getAtomPos(ring_atom_ids[i]) - center;
      RDGeom::Point3D v2 = conf.getAtomPos(ring_atom_ids[(i + 1) % ring_atom_ids.size()]) - center;
      norm += v1.crossProduct(v2);
    }
    norm /= static_cast<double>(ring_atom_ids.size());
    norm.normalize();
  }

  static void
  compute_normal_fast(const RDKit::RWMol *mol, const std::vector<int> &ring_, RDGeom::Point3D &norm1) {

    if (ring_.size() < 3) {
      throw std::invalid_argument("Ring must contain at least 3 atoms to compute a normal.");
    }

    const RDKit::Conformer &conf = mol->getConformer();

    std::vector<int> ring = ring_;
    std::sort(ring.begin(), ring.end());

    RDGeom::Point3D a = conf.getAtomPos(ring[0]);
    RDGeom::Point3D b = conf.getAtomPos(ring[1]);
    RDGeom::Point3D c = conf.getAtomPos(ring[2]);

    RDGeom::Point3D AB = b - a;
    RDGeom::Point3D AC = c - a;

    norm1 = AB.crossProduct(AC);
    norm1.normalize();
  }
};

} // namespace lahuta

#endif // LAHUTA_RINGS_HPP
