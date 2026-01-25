#ifndef LAHUTA_ENTITIES_RECORDS_HPP
#define LAHUTA_ENTITIES_RECORDS_HPP

//
// Previously, we were storing geometrical data here. This was changed when we enabled
// multi-conformer support to avoid stale geometry in these records. Now, geometry is
// computed on-the-fly via accessors that take a Conformer reference. - Besian, October 2025
//
#include <cmath>
#include <vector>

#include <Geometry/point.h>
#include <rdkit/GraphMol/Conformer.h>

#include "entities/entity_id.hpp"
#include "typing/types.hpp"

// clang-format off
namespace lahuta {

enum class FeatureGroup {
  None = 0,
  QuaternaryAmine,
  TertiaryAmine,
  Sulfonium,
  SulfonicAcid,
  Sulfate,
  Phosphate,
  Halocarbon,
  Guanidine,
  Acetamidine,
  Carboxylate
};

struct AtomRec {
  AtomType type;
  std::reference_wrapper<const RDKit::Atom>  atom; // Technically, atoms are stores as T*, but I'm not sure upstream will work
                                                   // with a nullable object. So we use T& as a mandatory never-null.
};

struct RingRec {
  std::vector<std::reference_wrapper<const RDKit::Atom>> atoms;
  bool aromatic;

  RDGeom::Point3D center(const RDKit::Conformer &conf) const {
    RDGeom::Point3D c{0.0, 0.0, 0.0};
    const auto n = atoms.size();
    if (n == 0) return c;
    for (const auto &a : atoms) c += conf.getAtomPos(a.get().getIdx());
    c /= static_cast<double>(n);
    if (n > 1) {
      // Clamp tiny values to zero to avoid floating-point artifacts from cancellation
      constexpr double CENTROID_EPSILON = 1e-6;
      if (std::abs(c.x) < CENTROID_EPSILON) c.x = 0.0;
      if (std::abs(c.y) < CENTROID_EPSILON) c.y = 0.0;
      if (std::abs(c.z) < CENTROID_EPSILON) c.z = 0.0;
    }
    return c;
  }

  RDGeom::Point3D normal(const RDKit::Conformer &conf) const {
    RDGeom::Point3D nrm{0.0, 0.0, 0.0};
    if (atoms.size() < 3) return nrm;
    const auto &a0 = atoms[0].get();
    const auto &a1 = atoms[1].get();
    const auto &a2 = atoms[2].get();
    auto p0 = conf.getAtomPos(a0.getIdx());
    auto p1 = conf.getAtomPos(a1.getIdx());
    auto p2 = conf.getAtomPos(a2.getIdx());
    auto v1 = p1 - p0;
    auto v2 = p2 - p0;
    nrm = v1.crossProduct(v2);
    const double len_sq = nrm.lengthSq();
    if (len_sq > 0.0) nrm /= std::sqrt(len_sq);
    return nrm;
  }
};

struct GroupRec {
  AtomType              a_type;
  FeatureGroup          type;
  std::vector<std::reference_wrapper<const RDKit::Atom>> atoms;

  RDGeom::Point3D center(const RDKit::Conformer &conf) const {
    RDGeom::Point3D c{0.0, 0.0, 0.0};
    const auto n = atoms.size();
    if (n == 0) return c;
    for (const auto &a : atoms) c += conf.getAtomPos(a.get().getIdx());
    c /= static_cast<double>(n);
    if (n > 1) {
      // Clamp tiny values to zero to avoid floating-point artifacts from cancellation
      constexpr double CENTROID_EPSILON = 1e-6;
      if (std::abs(c.x) < CENTROID_EPSILON) c.x = 0.0;
      if (std::abs(c.y) < CENTROID_EPSILON) c.y = 0.0;
      if (std::abs(c.z) < CENTROID_EPSILON) c.z = 0.0;
    }
    return c;
  }
};

template<typename T> struct KindOf;
template<> struct KindOf<AtomRec>  { static constexpr Kind value = Kind::Atom; };
template<> struct KindOf<RingRec>  { static constexpr Kind value = Kind::Ring; };
template<> struct KindOf<GroupRec> { static constexpr Kind value = Kind::Group; };

template<Kind K> struct RecordTypeFor;
template<> struct RecordTypeFor<Kind::Atom>  { using type = AtomRec; };
template<> struct RecordTypeFor<Kind::Ring>  { using type = RingRec; };
template<> struct RecordTypeFor<Kind::Group> { using type = GroupRec; };

} // namespace lahuta

#endif // LAHUTA_ENTITIES_RECORDS_HPP
