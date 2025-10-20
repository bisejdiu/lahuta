#ifndef LAHUTA_RINGS_HPP
#define LAHUTA_RINGS_HPP

#include <vector>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/RWMol.h>

namespace lahuta {

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

  static void compute_normal(const RDKit::RWMol *mol, const std::vector<int> &ring_atom_ids, const RDGeom::Point3D &center, RDGeom::Point3D &norm) {

    if (ring_atom_ids.size() < 3) {
      throw std::invalid_argument("Ring must contain at least 3 atoms to compute a normal.");
    }

    const RDKit::Conformer &conf = mol->getConformer();

    norm = RDGeom::Point3D{0.0, 0.0, 0.0};
    for (size_t i = 0; i < ring_atom_ids.size(); ++i) {
      RDGeom::Point3D v1 = conf.getAtomPos(ring_atom_ids[i]) - center;
      RDGeom::Point3D v2 = conf.getAtomPos(ring_atom_ids[(i + 1) % ring_atom_ids.size()]) - center;

      auto new_norm = v1.crossProduct(v2);
      new_norm.normalize();

      norm += new_norm;
    }
    norm /= static_cast<double>(ring_atom_ids.size());
  }

  static void compute_normal_fast(const RDKit::RWMol *mol, const std::vector<int> &ring_, RDGeom::Point3D &norm1) {

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
