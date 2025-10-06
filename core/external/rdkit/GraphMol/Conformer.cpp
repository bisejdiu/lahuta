//
//  Copyright (C) 2001-2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Conformer.h"
#include "ROMol.h"

namespace RDKit {

void Conformer::setOwningMol(ROMol *mol) {
  PRECONDITION(mol, "");
  dp_mol = mol;
}

void Conformer::setOwningMol(ROMol &mol) { setOwningMol(&mol); }

const RDGeom::POINT3D_VECT &Conformer::getPositions() const {
  const auto &pos = d_external_positions ? *d_external_positions : d_positions;
  if (dp_mol) {
    PRECONDITION(dp_mol->getNumAtoms() == pos.size(), "");
  }
  return pos;
}

RDGeom::POINT3D_VECT &Conformer::getPositions() {
  if (d_external_positions) {
    throw ConformerException("getPositions(): Mutation is not allowed when external coordinate view is bound");
  }
  return d_positions;
}

const RDGeom::Point3D &Conformer::getAtomPos(unsigned int atomId) const {
  const auto &pos = d_external_positions ? *d_external_positions : d_positions;
  if (dp_mol) {
    PRECONDITION(dp_mol->getNumAtoms() == pos.size(), "");
  }
  URANGE_CHECK(atomId, pos.size());
  return pos.at(atomId);
}

RDGeom::Point3D &Conformer::getAtomPos(unsigned int atomId) {
  if (d_external_positions) {
    throw ConformerException("getAtomPos(): Mutation is not allowed when external coordinate view is bound");
  }
  if (dp_mol) {
    PRECONDITION(dp_mol->getNumAtoms() == d_positions.size(), "");
  }
  URANGE_CHECK(atomId, d_positions.size());
  return d_positions.at(atomId);
}
}  // namespace RDKit
