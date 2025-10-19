#ifndef LAHUTA_BOND_ORDER_HPP
#define LAHUTA_BOND_ORDER_HPP

#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>

// clang-format off
using namespace RDKit;
namespace lahuta {

using HybridizationType = RDKit::Atom::HybridizationType;
using SubStrMatches = std::vector<RDKit::MatchVectType>;

constexpr std::pair<const char *, HybridizationType> smartsList[] = {
    {"[D4]", HybridizationType::SP3},
    {"[D5]", HybridizationType::SP3D},
    {"[D6]", HybridizationType::SP3D2},
    {"[C]", HybridizationType::SP3},
    {"[c,$(C=*)]", HybridizationType::SP2},
    {"[$(C#*),$(C(=*)=*)]", HybridizationType::SP},

    {"[N]", HybridizationType::SP3},
    {"[n,$(N=*),$(N[#6,#7,#8]=,:,#*)]", HybridizationType::SP2},
    {"[ND1,ND2,ND3]a", HybridizationType::SP2},
    {"[$(N#*),$([ND2](=*)=*)]", HybridizationType::SP},

    {"[O]", HybridizationType::SP3},
    {"[o,$(O=*),$(O[#6,#7,#8]=,:*)]", HybridizationType::SP2},
    {"[$([#8D1][#6][#8D1])]", HybridizationType::SP2},
    {"[$(O#*)]", HybridizationType::SP},

    {"[S]", RDKit::Atom::HybridizationType::SP3},
    {"[#16;s,$([SD1]=*)]", RDKit::Atom::HybridizationType::SP2},
    {"[SD6]", RDKit::Atom::HybridizationType::SP3D2},
};

using SmartsPatternPair = std::pair<const char *, std::vector<int>>;
const std::vector<SmartsPatternPair> bondSmarts = {
    {"[x2,x3]1[#6]([#7D3]2)[#6][#6][#6]2[x2,x3][#6]([#7D3]3)[#6][#6][#6]3["
     "x2,"
     "x3][#6]([#7D3]4)[#6][#6][#6]4[x2,x3][#6]([#7D3]5)[#6][#6][#6]51",
     {0,  1,  2,  1,  2,  1,  1,  3,  1,  3,  4,  2,  4,  5,  1,  5,  2,
      1,  5,  6,  2,  6,  7,  1,  7,  8,  2,  7,  9,  1,  9,  10, 2,  10,
      11, 1,  11, 8,  1,  11, 12, 2,  12, 13, 1,  13, 14, 1,  13, 15, 2,
      15, 16, 1,  16, 17, 2,  17, 14, 1,  17, 18, 1,  18, 19, 2,  19, 20,
      1,  19, 21, 1,  21, 22, 2,  22, 23, 1,  23, 20, 2}},
    {"[x2,x3]1[#6]([#7D3]2)[#6][#6][#6]2[x2,x3][#6]([#7]3)[#6][#6][#6]3[x2,"
     "x3][#6]([#7D3]4)[#6][#6][#6]4[x2,x3][#6]([#7]5)[#6][#6][#6]51",
     {0,  1,  2,  1,  2,  1,  1,  3,  1,  3,  4,  2,  4,  5,  1,  5,  2,
      1,  5,  6,  2,  6,  7,  1,  7,  8,  2,  7,  9,  1,  9,  10, 2,  10,
      11, 1,  11, 8,  1,  11, 12, 2,  12, 13, 1,  13, 14, 1,  13, 15, 2,
      15, 16, 1,  16, 17, 2,  17, 14, 1,  17, 18, 1,  18, 19, 2,  19, 20,
      1,  19, 21, 1,  21, 22, 2,  22, 23, 1,  23, 20, 2}},
    {"[x2,x3]1[#6]([#7]2)[#6][#6][#6]2[x2,x3][#6]([#7]3)[#6][#6][#6]3[x2,"
     "x3]["
     "#6]([#7]4)[#6][#6][#6]4[x2,x3][#6]([#7]5)[#6][#6][#6]51",
     {0,  1,  2,  1,  2,  1,  1,  3,  1,  3,  4,  2,  4,  5,  1,  5,  2,
      1,  5,  6,  2,  6,  7,  1,  7,  8,  2,  7,  9,  1,  9,  10, 2,  10,
      11, 1,  11, 8,  1,  11, 12, 2,  12, 13, 1,  13, 14, 1,  13, 15, 2,
      15, 16, 1,  16, 17, 2,  17, 14, 1,  17, 18, 1,  18, 19, 2,  19, 20,
      1,  19, 21, 1,  21, 22, 2,  22, 23, 1,  23, 20, 2}},
    {"[#7D2][#7D2^1][#7D1]", {0, 1, 2, 1, 2, 2}},
    {"[#8D1][#7D3^2]([#8D1])*", {0, 1, 2, 1, 2, 2, 1, 3, 1}},
    {"[#16D4]([#8D1])([#8D1])([*!#8])([*!#8])",
     {0, 1, 2, 0, 2, 2, 0, 3, 1, 0, 4, 1}},
    {"[#16D4]([#8D1])([#8D1])([#8-,#8D1])([#8-,#8D1])",
     {0, 1, 2, 0, 2, 2, 0, 3, 1, 0, 4, 1}},
    {"[#16D4]([#16D1])([#8D1])([#8-,#8])([#8-,#8])",
     {0, 1, 2, 0, 2, 2, 0, 3, 1, 0, 4, 1}},
    {"[#16D3]([#8D1])([*!#8])([*!#8])", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#16D3]([#8D1])([#8D1-])([#8D1-])", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#16D3]([#8D1])([#8])([#8])", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#16D2]([#8D1])([#16D1])", {0, 1, 2, 0, 2, 2}},
    {"[#16D2]([#8D1])([*!#8])", {0, 1, 2, 0, 2, 1}},
    {"[#16D2]([#8D1])([#8D1])", {0, 1, 2, 0, 2, 2}},
    {"[#15D3]([#8D1])([#8D1])([#8D2])", {0, 1, 2, 0, 2, 2, 0, 3, 1}},
    {"[#7D2]([#8D1])([#1])", {0, 1, 2, 0, 2, 1}},
    {"[#15D4]([#8D1])(*)(*)(*)", {0, 1, 2, 0, 2, 1, 0, 3, 1, 0, 4, 1}},
    {"[#6D3^2]([#8D1])([#8])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#8D1][#6D2^1][#8D1]", {0, 1, 2, 1, 2, 2}},
    {"[#6D3^2]([#8D1;!-])([#7])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#34D3^2]([#8D1])([#8])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#6D3^2]([#8D1])([#16])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#6D3^2]([#16D1])([#16])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[CD3^2]([#16D1])([N])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
};

inline void OBBondTypeAssignment(RDKit::ROMol &mol) {

  // Precompute patterns statically
  static std::vector<RDKit::ROMol *> patterns = [] {
    std::vector<RDKit::ROMol *> temp;
    temp.resize(bondSmarts.size());
    for (size_t i = 0; i < bondSmarts.size(); ++i) {
      temp[i] = RDKit::SmartsToMol(bondSmarts[i].first);
    }
    return temp;
  }();

  for (size_t i = 0; i < bondSmarts.size(); ++i) {
    const auto &[smarts, bondVector] = bondSmarts[i];
    RDKit::ROMol *pattern = patterns[i];

    SubstructMatchParameters params;
    params.maxMatches = 100000;
    std::vector<RDKit::MatchVectType> matchList;
    matchList = RDKit::SubstructMatch(mol, *pattern, params);

    for (const auto &match : matchList) {
      for (auto j = 0; j < bondVector.size(); j += 3) {
        auto bond = mol.getBondBetweenAtoms(match[bondVector[j]].second, match[bondVector[j + 1]].second);
        if (bond) {
          bond->setBondType(static_cast<RDKit::Bond::BondType>(bondVector[j + 2]));
        }
      }
    }
  }
}

inline double average_bond_angle(const RDKit::Atom *atom) {
  double avgDegrees = 0.0;
  int n = 0;

  const RDKit::ROMol &mol = atom->getOwningMol();
  const RDKit::Conformer &conf = mol.getConformer();

  for (auto bondIt1 = mol.getAtomBonds(atom); bondIt1.first != bondIt1.second;
       ++bondIt1.first) {
    const RDKit::Bond *bond1 = mol[*bondIt1.first];
    const RDKit::Atom *neighbor1 = bond1->getOtherAtom(atom);

    for (auto bondIt2 = bondIt1; ++bondIt2.first != bondIt2.second;) {
      const RDKit::Bond *bond2 = mol[*bondIt2.first];
      const RDKit::Atom *neighbor2 = bond2->getOtherAtom(atom);

      // Get the coordinates of the atoms
      RDGeom::Point3D pos1 = conf.getAtomPos(atom->getIdx());
      RDGeom::Point3D pos2 = conf.getAtomPos(neighbor1->getIdx());
      RDGeom::Point3D pos3 = conf.getAtomPos(neighbor2->getIdx());

      // Calculate the vectors
      RDGeom::Point3D v1 = pos2 - pos1;
      RDGeom::Point3D v2 = pos3 - pos1;

      // Normalize the vectors
      v1.normalize();
      v2.normalize();

      // Calculate the angle in radians
      double angle = acos(v1.dotProduct(v2));

      // Convert to degrees
      double degrees = angle * 180.0 / M_PI;

      avgDegrees += degrees;
      n++;
    }
  }

  if (n >= 1) {
    avgDegrees /= n;
  }

  return avgDegrees;
}

/*!  This function calculates the torsion angle of three vectors, represented
  by four points A--B--C--D, i.e. B and C are vertexes, but none of A--B,
  B--C, and C--D are colinear.  A "torsion angle" is the amount of "twist"
  or torsion needed around the B--C axis to bring A--B into the same plane
  as B--C--D.  The torsion is measured by "looking down" the vector B--C so
  that B is superimposed on C, then noting how far you'd have to rotate
  A--B to superimpose A over D.  Angles are + in the anti-clockwise
  direction.  The operation is symmetrical in that if you reverse the image
  (look from C to B and rotate D over A), you get the same answer.
*/
inline double compute_dihedral(const RDGeom::Point3D &a, const RDGeom::Point3D &b,
                        const RDGeom::Point3D &c, const RDGeom::Point3D &d) {

  double torsion;
  RDGeom::Point3D b1, b2, b3, c1, c2, c3;

  b1 = a - b;
  b2 = b - c;
  b3 = c - d;

  double rb2 = sqrt(b2.dotProduct(b2));

  RDGeom::Point3D b2xb3 = b2.crossProduct(b3);
  RDGeom::Point3D b1xb2 = b1.crossProduct(b2);
  torsion = -atan2(rb2 * b1.dotProduct(b2xb3), b1xb2.dotProduct(b2xb3));

  return (torsion * 180.0 / M_PI);
}


void perceive_bond_orders_obabel(RDKit::RWMol &mol);

} // namespace lahuta

#endif // LAHUTA_BOND_ORDER_HPP
