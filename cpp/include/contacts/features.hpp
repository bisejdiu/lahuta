#ifndef LAHUTA_FEATURES_HPP
#define LAHUTA_FEATURES_HPP

#include "atom_types.hpp"

namespace lahuta {

class Luni;
using FeatureTypeCheckFunc = std::function<bool(const AtomType &, const AtomType &)>;

enum class FeatureGroup {
  None = 0,
  QuaternaryAmine = 1,
  TertiaryAmine = 2,
  Sulfonium = 3,
  SulfonicAcid = 4,
  Sulfate = 5,
  Phosphate = 6,
  Halocarbon = 7,
  Guanidine = 8,
  Acetamidine = 9,
  Carboxylate = 10
};

class FeatureVec;

struct Feature {
  AtomType type;
  FeatureGroup group;
  std::vector<const RDKit::Atom *> members;
  RDGeom::Point3D center;

  int get_id() const { return id; }

  Feature() = default;
  Feature(AtomType type, FeatureGroup group, std::vector<const RDKit::Atom *> members, RDGeom::Point3D center)
      : type(type), group(group), members(std::move(members)), center(center) {}

private:
  Feature(
      AtomType type, FeatureGroup group, std::vector<const RDKit::Atom *> members, RDGeom::Point3D center,
      int id)
      : type(type), group(group), members(std::move(members)), center(center), id(id) {}

public:
  // only GroupTypeStrategy can set the id
  friend class GroupTypeStrategy;
  friend FeatureVec get_features(const Luni *luni, AtomType type, FeatureTypeCheckFunc check_func);

private:
  int id;
};

struct FeatureVec {
  std::vector<Feature> features;

  Feature &operator[](size_t index) { return features[index]; }

  const Feature &operator[](size_t index) const { return features[index]; }

  const RDGeom::POINT3D_VECT *centers() const {
    auto *center_vec = new RDGeom::POINT3D_VECT();
    center_vec->reserve(features.size());
    for (const auto &feature : features) {
      center_vec->push_back(feature.center);
    }
    return center_vec;
  }

  std::vector<const RDKit::Atom *> atoms() const {
    std::vector<const RDKit::Atom *> atoms_vec;
    atoms_vec.reserve(features.size());
    for (const auto &feature : features) {
      atoms_vec.insert(atoms_vec.end(), feature.members.begin(), feature.members.end());
    }
    return atoms_vec;
  }

  RDGeom::POINT3D_VECT positions() const {
    if (features.empty()) {
      return {};
    }
    RDGeom::POINT3D_VECT pos_vec;
    const auto &conf = features.front().members.front()->getOwningMol().getConformer();
    for (const auto &feature : features) {
      pos_vec.push_back(feature.center);
    }
    return pos_vec;
  }
};

/*using FeatureGroupCheckFunc = std::function<bool(const AtomType &, const FeatureGroup &)>;*/
FeatureVec
get_features(const Luni *luni, AtomType type, FeatureTypeCheckFunc check_func = AtomTypeFlags::has_any);

} // namespace lahuta

#endif // LAHUTA_FEATURES_HPP
