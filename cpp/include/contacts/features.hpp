#ifndef LAHUTA_FEATURES_HPP
#define LAHUTA_FEATURES_HPP

#include "atom_types.hpp"

namespace lahuta {

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

struct Feature {
  AtomType type;
  FeatureGroup group;
  std::vector<const RDKit::Atom *> members;
  double center[3];

  int get_id() const { return id; }

  // only GroupTypeStrategy can set the id
  friend class GroupTypeStrategy;

private:
  int id;
};

} // namespace lahuta

#endif // LAHUTA_FEATURES_HPP
