//
//  Copyright (C) 2001-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_SMILESPARSEOPS_H
#define RD_SMILESPARSEOPS_H
#include <GraphMol/Bond.h>

namespace RDKit {
class RWMol;
class Atom;
class QueryBond;
}  // namespace RDKit
namespace SmilesParseOps {
void CheckRingClosureBranchStatus(RDKit::Atom *atom,
                                                           RDKit::RWMol *mp);
void ReportParseError(const char *message,
                                               bool throwIt = true);
void CleanupAfterParseError(RDKit::RWMol *mol);
// This uses SMARTS semantics: unspecified bonds are treated as
// aromatic or single.
void AddFragToMol(
    RDKit::RWMol *mol, RDKit::RWMol *frag,
    RDKit::Bond::BondType bondOrder = RDKit::Bond::UNSPECIFIED,
    RDKit::Bond::BondDir bondDir = RDKit::Bond::NONE);
RDKit::Bond::BondType GetUnspecifiedBondType(
    const RDKit::RWMol *mol, const RDKit::Atom *atom1,
    const RDKit::Atom *atom2);
void CheckChiralitySpecifications(RDKit::RWMol *mol,
                                                           bool strict);
void CloseMolRings(RDKit::RWMol *mol,
                                            bool toleratePartials);
void SetUnspecifiedBondTypes(RDKit::RWMol *mol);
void AdjustAtomChiralityFlags(RDKit::RWMol *mol);
void CleanupAfterParsing(RDKit::RWMol *mol);
void parseCXExtensions(
    RDKit::RWMol &mol, const std::string &extText,
    std::string::const_iterator &pos, unsigned int startAtomIdx = 0,
    unsigned int startBondIdx = 0);
inline void parseCXExtensions(RDKit::RWMol &mol, const std::string &extText,
                              unsigned int startAtomIdx,
                              unsigned int startBondIdx) {
  auto iter = extText.begin();
  parseCXExtensions(mol, extText, iter, startAtomIdx, startBondIdx);
};
//! removes formal charge, isotope, etc. Primarily useful for QueryAtoms
void ClearAtomChemicalProps(RDKit::Atom *atom);

//! returns whether or not the combination of tag and permutation provided are
//! legal
bool checkChiralPermutation(int chiralTag,
                                                     int permutation);

//! this is a bit of a hack to try and get nicer "SMILES" from
//! a SMARTS molecule
RDKit::QueryBond *getUnspecifiedQueryBond(
    const RDKit::Atom *a1, const RDKit::Atom *a2);

namespace detail {
constexpr auto _needsDetectBondStereo = "_needsDetectBondStereo";
constexpr auto _needsDetectAtomStereo = "_needsDetectAtomStereo";
}  // namespace detail
}  // namespace SmilesParseOps

#endif
