#include <GraphMol/Atom.h>
#include <GraphMol/RWMol.h>

inline unsigned int ob_explicit_valence(const RDKit::RWMol &mol,
                                             const RDKit::Atom *atom) {
  unsigned int valence = 0;
  for (auto bondIt = mol.getAtomBonds(atom); bondIt.first != bondIt.second;
      ++bondIt.first) {
      valence++;
  }
  return valence;
}

inline double bond_length_sq(const RDKit::Conformer &conf,
                             const RDKit::Bond *bond) {
  auto a1Pos = conf.getAtomPos(bond->getBeginAtomIdx());
  auto a2Pos = conf.getAtomPos(bond->getEndAtomIdx());

  double d2x = (a1Pos.x - a2Pos.x);
  double d2y = (a1Pos.y - a2Pos.y);
  double d2z = (a1Pos.z - a2Pos.z);

  return d2x * d2x + d2y * d2y + d2z * d2z;
};

inline double compute_angle(const RDKit::Conformer &conf,
                            const RDKit::Atom *atom1, const RDKit::Atom *atom2,
                            RDKit::Atom *atom3) {

  RDGeom::Point3D v1, v2;
  RDGeom::Point3D p1 = conf.getAtomPos(atom1->getIdx());
  RDGeom::Point3D p2 = conf.getAtomPos(atom2->getIdx());
  RDGeom::Point3D p3 = conf.getAtomPos(atom3->getIdx());

  v1 = p1 - p2;
  v2 = p3 - p2;

  if (v1.length() < 0.0001 || v2.length() < 0.0001) {
    return 0.0;
  }

  double angle = v1.angleTo(v2) * 180.0 / M_PI;

  return angle;
}

inline double smallest_bond_angle(const RDKit::RWMol &mol,
                                  const RDKit::Conformer &conf,
                                  const RDKit::Atom *atom) {
  double min_degrees = 360.0;

  for (auto it = mol.getAtomBonds(atom); it.first != it.second; ++it.first) {
    auto bond1 = mol[*it.first];
    auto atom1 = bond1->getOtherAtom(atom);
    for (auto it2 = mol.getAtomBonds(atom); it2.first != it2.second;
         ++it2.first) {
      if (it.first == it2.first)
        continue;
      auto bond2 = mol[*it2.first];
      auto atom2 = bond2->getOtherAtom(atom);

      double angle = compute_angle(conf, atom1, atom, atom2);
      if (angle < min_degrees) {
        min_degrees = angle;
      }
    }
  }

  return min_degrees;
}

void clean_bonds(RDKit::RWMol &mol, RDKit::Conformer &conf);
