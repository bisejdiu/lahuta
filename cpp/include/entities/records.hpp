#ifndef LAHUTA_ENTITIES_RECORDS_HPP
#define LAHUTA_ENTITIES_RECORDS_HPP

#include "Geometry/point.h"
#include "typing/types.hpp"
#include <entities/entity_id.hpp>
#include <vector>

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
  const RDKit::Atom &atom; // Technically, atoms are stores as T*, but I'm not sure upstream will work
                           // with a nullable object. So we use T& as a mandatory never-null.
};

struct RingRec {
  std::vector<std::reference_wrapper<const RDKit::Atom>> atoms;
  RDGeom::Point3D       center;
  RDGeom::Point3D       normal;
  bool aromatic;
};

struct GroupRec {
  AtomType              a_type;
  FeatureGroup          type;
  // std::vector<uint32_t> atoms;
  std::vector<std::reference_wrapper<const RDKit::Atom>> atoms;
  RDGeom::Point3D       center;
};

template <typename T> struct KindOf;
template <> struct KindOf<AtomRec>  { static constexpr Kind value = Kind::Atom; };
template <> struct KindOf<RingRec>  { static constexpr Kind value = Kind::Ring; };
template <> struct KindOf<GroupRec> { static constexpr Kind value = Kind::Group; };

template <Kind K> struct RecordTypeFor;
template <> struct RecordTypeFor<Kind::Atom>  { using type = AtomRec; };
template <> struct RecordTypeFor<Kind::Ring>  { using type = RingRec; };
template <> struct RecordTypeFor<Kind::Group> { using type = GroupRec; };

} // namespace lahuta

#endif // LAHUTA_ENTITIES_RECORDS_HPP
