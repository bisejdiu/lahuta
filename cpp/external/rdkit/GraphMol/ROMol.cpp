//
//  Copyright (C) 2003-2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <boost/graph/compressed_sparse_row_graph.hpp>

// our stuff
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include "ROMol.h"
#include "Atom.h"
#include "Bond.h"
#include "Conformer.h"
#include "SubstanceGroup.h"
#include "AtomIterators.h"
#include "BondIterators.h"

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <RDGeneral/BoostEndInclude.h>
#endif

namespace RDKit {
class QueryAtom;
class QueryBond;

ROMol::AtomIterator ROMol::beginAtoms() { return AtomIterator(this); }
ROMol::ConstAtomIterator ROMol::beginAtoms() const {
  return ConstAtomIterator(this);
}
ROMol::AtomIterator ROMol::endAtoms() {
  return AtomIterator(this, getNumAtoms());
}
ROMol::ConstAtomIterator ROMol::endAtoms() const {
  return ConstAtomIterator(this, getNumAtoms());
}

ROMol::BondIterator ROMol::beginBonds() { return BondIterator(this); }
ROMol::ConstBondIterator ROMol::beginBonds() const {
  return ConstBondIterator(this);
}
ROMol::BondIterator ROMol::endBonds() {
  auto [beg, end] = getEdges();
  return BondIterator(this, end);
}
ROMol::ConstBondIterator ROMol::endBonds() const {
  auto [beg, end] = getEdges();
  return ConstBondIterator(this, end);
}

#ifdef RDK_USE_BOOST_SERIALIZATION
template <class Archive>
void ROMol::save(Archive &ar, const unsigned int) const {
  std::string pkl;
  MolPickler::pickleMol(*this, pkl, PicklerOps::AllProps);
  ar << pkl;
}

template <class Archive>
void ROMol::load(Archive &ar, const unsigned int) {
  std::string pkl;
  ar >> pkl;

  delete dp_ringInfo;
  initMol();

  numBonds = 0;
  MolPickler::molFromPickle(pkl, *this, PicklerOps::AllProps);
  numBonds = rdcast<unsigned int>(boost::num_edges(d_graph));
}

template RDKIT_GRAPHMOL_EXPORT void ROMol::save<boost::archive::text_oarchive>(
    boost::archive::text_oarchive &, const unsigned int) const;
template RDKIT_GRAPHMOL_EXPORT void ROMol::load<boost::archive::text_iarchive>(
    boost::archive::text_iarchive &, const unsigned int);
#endif

}  // namespace RDKit
