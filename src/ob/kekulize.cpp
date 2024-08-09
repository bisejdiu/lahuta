/**********************************************************************
kekulize.cpp - Kekulize a molecule

Copyright (C) 2017 Noel M. O'Boyle

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
/*
 * File:   kekulize.cpp
 * Modified by:  Besian I. Sejdiu
 * Redistributed as part of the Lahuta project and in compliance with the
 * GNU General Public License. See Open Babel (https://openbabel.org/)
 * for more information. The Open Babel license information is included
 * with the redistributed code ($LAHUTA_ROOT/exteranl/ob/LICENSE) and the
 * Lahuta license information is $LAHUTA_ROOT/LICENSE.
 */
// #include <boost/range/iterator_range_core.hpp>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Bond.h>
#include <rdkit/GraphMol/RWMol.h>

#include <cstdlib>
#include <cstring>
#include "ob/bitvec.h"
#include "ob/kekulize.h"

using namespace RDKit;
// using namespace OB;

static unsigned int GetMaxAtomIdx(RWMol *mol) { return mol->getNumAtoms() + 1; }
static unsigned int GetMaxBondIdx(RWMol *mol) { return mol->getNumBonds(); }

class Kekulizer {
public:
  Kekulizer(RWMol *mol)
      : m_mol(mol), needs_dbl_bond(nullptr), doubleBonds(nullptr),
        kekule_system(nullptr) {
    atomArraySize = GetMaxAtomIdx(m_mol) + 1;
    bondArraySize = GetMaxBondIdx(m_mol) + 1;
  }
  ~Kekulizer() {
    delete needs_dbl_bond;
    delete doubleBonds;
    delete kekule_system;
  }
  bool GreedyMatch();
  bool BackTrack();
  void AssignDoubleBonds();

private:
  bool FindPath(unsigned int atomidx, bool isDoubleBond, OBBitVec &visited);
  RWMol *m_mol;
  OBBitVec *needs_dbl_bond;
  OBBitVec *doubleBonds;
  OBBitVec *kekule_system;
  unsigned int atomArraySize;
  unsigned int bondArraySize;
  std::vector<unsigned int> m_path;
};

static bool IsSpecialCase(Atom *atom) {
  switch (atom->getAtomicNum()) {
  case 7:
    // Any exo-cyclic double bond from a N
    // e.g. pyridine N-oxide as the double bond form
    if (atom->getTotalDegree() == 3 && atom->getFormalCharge() == 0)
      return true;
    break;
  case 16: // e.g. Cs1(=O)ccccn1 but not O=s1(=O)cccn1
    if (atom->getTotalDegree() == 4 && atom->getFormalCharge() == 0 &&
        atom->getTotalValence() < 6)
      return true;
    break;
  }
  return false;
}

static bool NeedsDoubleBond(Atom *atom) {
  if (!atom->getIsAromatic())
    return false;

  // Does it already have an explicit double bond?
  for (const auto &bondItr :
       boost::make_iterator_range(atom->getOwningMol().getAtomBonds(atom))) {
    const Bond *bond = atom->getOwningMol()[bondItr];
    if (bond->getIsAromatic())
      continue;
    const Atom *nbr = bond->getOtherAtom(atom);
    switch (bond->getBondType()) {
    case Bond::UNSPECIFIED:
    case Bond::SINGLE:
      continue;
    case Bond::DOUBLE:
      if (IsSpecialCase(atom))
        return true;
      return false;
    default: // bond order > 2
      return false;
    }
  }

  // Is it one of the cases where we know that it only has single bonds?
  // NOTE: For RDKit we need to use getDegree() instead of getTotalDegree()
  int chg = atom->getFormalCharge();
  int deg = atom->getDegree();
  switch (atom->getAtomicNum()) {
  case 6:
    if (deg == 3 && (chg == 1 || chg == -1))
      return false;
    break;
  case 5:
  case 7:
  case 15:
  case 33:
  case 51:
  case 83:
    switch (chg) {
    case 0: // e.g. a pyrrole-type nitrogen
      if (deg == 3 || deg > 4)
        return false;
      break;
    case -1:
      if (deg == 2)
        return false;
      break;
    case 1:
      if (deg > 3)
        return false;
    }
    break;
  case 8:
  case 16:
  case 34:
  case 52:
    switch (chg) {
    case 0:
      if (deg == 2 || deg == 4 || deg > 5)
        return false;
      break;
    case -1:
    case 1:
      if (deg == 3 || deg == 5 || deg > 6)
        return false;
    }
  }

  return true; // It needs a double bond
}

class NodeIterator {
public:
  NodeIterator(unsigned int *&degrees, unsigned int atomArraySize)
      : m_degrees(degrees), m_atomArraySize(atomArraySize), m_counter(0),
        finishedDegTwo(false) {}
  unsigned int next() {
    m_counter++;

    if (!finishedDegTwo) { // return deg 2 nodes first
      for (; m_counter < m_atomArraySize; ++m_counter) {
        if (m_degrees[m_counter] == 2) {
          return m_counter;
        }
      }
      finishedDegTwo = true;
      m_counter = 1; // first atom has idx 1
    }

    // return nodes with degree > 2
    for (; m_counter < m_atomArraySize; ++m_counter) {
      if (m_degrees[m_counter] > 2) {
        return m_counter;
      }
    }

    // Finished - return 0 signalling the end of iteration
    return 0;
  }

private:
  unsigned int *&m_degrees;
  unsigned int m_atomArraySize;
  unsigned int m_counter;
  bool finishedDegTwo;
};

void Kekulizer::AssignDoubleBonds() {
  int bit;
  for (bit = doubleBonds->FirstBit(); bit != doubleBonds->EndBit();
       bit = doubleBonds->NextBit(bit)) {
    m_mol->getBondWithIdx(bit)->setBondType(Bond::DOUBLE);
  }
}

bool Kekulizer::GreedyMatch() {

  // What atoms need a double bond? The job of kekulization is
  // to give all of these atoms a single double bond.
  needs_dbl_bond = new OBBitVec(atomArraySize); // defaults to all False
  for (auto *atom : m_mol->atoms()) {
    if (NeedsDoubleBond(atom)) {
      needs_dbl_bond->SetBitOn(atom->getIdx());
    }
  }
  // Make a copy of needs_dbl_bond, to restrict the traversal in BackTrack()
  kekule_system = new OBBitVec(*needs_dbl_bond);

  // Create lookup of degrees
  unsigned int *degrees =
      (unsigned int *)calloc(atomArraySize, sizeof(unsigned int));
  std::vector<Atom *> degreeOneAtoms;
  for (auto *atom : m_mol->atoms()) {
    unsigned int atom_idx = atom->getIdx();
    if (!needs_dbl_bond->BitIsSet(atom_idx)) {
      degrees[atom_idx] = 0;
      continue;
    }
    unsigned int mdeg = 0;
    for (const auto &bondItr :
         boost::make_iterator_range(atom->getOwningMol().getAtomBonds(atom))) {
      const Bond *bond = atom->getOwningMol()[bondItr];
      if (!bond->getIsAromatic())
        continue;
      const Atom *nbr = bond->getOtherAtom(atom);
      if (needs_dbl_bond->BitIsSet(nbr->getIdx()))
        mdeg++;
    }
    degrees[atom_idx] = mdeg;
    if (mdeg == 1)
      degreeOneAtoms.push_back(atom);
  }

  // Location of assigned double bonds
  doubleBonds = new OBBitVec(bondArraySize); // defaults to all False

  bool finished = false;
  while (true) { // Main loop

    // Complete all of the degree one nodes
    while (!degreeOneAtoms.empty()) {
      Atom *atom = degreeOneAtoms.back();
      degreeOneAtoms.pop_back();
      // some nodes may already have been handled
      if (!needs_dbl_bond->BitIsSet(atom->getIdx()))
        continue;
      for (const auto &bondItr : boost::make_iterator_range(
               atom->getOwningMol().getAtomBonds(atom))) {
        const Bond *bond = atom->getOwningMol()[bondItr];
        if (!bond->getIsAromatic())
          continue;
        Atom *nbr = bond->getOtherAtom(atom);
        if (!needs_dbl_bond->BitIsSet(nbr->getIdx()))
          continue;
        doubleBonds->SetBitOn(bond->getIdx());
        needs_dbl_bond->SetBitOff(atom->getIdx());
        needs_dbl_bond->SetBitOff(nbr->getIdx());
        // now update degree information for nbr's neighbors
        for (const auto &nbrbondItr : boost::make_iterator_range(
                 nbr->getOwningMol().getAtomBonds(nbr))) {
          const Bond *nbrbond = nbr->getOwningMol()[nbrbondItr];
          if (&*(nbrbond->getOtherAtom(nbr)) == &*atom ||
              !nbrbond->getIsAromatic())
            continue;
          Atom *nbrnbr = nbrbond->getOtherAtom(nbr);
          unsigned int nbrnbrIdx = nbrnbr->getIdx();
          if (!needs_dbl_bond->BitIsSet(nbrnbrIdx))
            continue;
          degrees[nbrnbrIdx]--;
          if (degrees[nbrnbrIdx] == 1)
            degreeOneAtoms.push_back(nbrnbr);
        }
        break;
      }
    }

    if (needs_dbl_bond->IsEmpty()) {
      finished = true;
      break;
    }

    // Now handle any remaining degree 2 or 3 nodes
    // We handle deg 2 nodes first and then 3, and the iteration over these
    // nodes is abstracted away. Once a double-bond is added that generates more
    // degree one nodes, then the iterator is exited
    NodeIterator iterator(degrees, atomArraySize);
    bool change = false;
    while (unsigned int atomIdx = iterator.next()) {
      if (!needs_dbl_bond->BitIsSet(atomIdx))
        continue;
      // The following is almost identical to the code above for deg 1 atoms
      // except for handling the variable 'change'
      Atom *atom = m_mol->getAtomWithIdx(atomIdx);
      for (const auto &bondItr : boost::make_iterator_range(
               atom->getOwningMol().getAtomBonds(atom))) {
        Bond *bond = atom->getOwningMol()[bondItr];
        if (!bond->getIsAromatic())
          continue;
        Atom *nbr = bond->getOtherAtom(atom);
        if (!needs_dbl_bond->BitIsSet(nbr->getIdx()))
          continue;
        // create a double bond from atom -> nbr
        doubleBonds->SetBitOn(bond->getIdx());
        needs_dbl_bond->SetBitOff(atomIdx);
        needs_dbl_bond->SetBitOff(nbr->getIdx());
        // now update degree information for both atom's and nbr's neighbors
        for (int N = 0; N < 2; N++) {
          Atom *ref = N == 0 ? atom : nbr;
          for (const auto &nbrbondItr : boost::make_iterator_range(
                   ref->getOwningMol().getAtomBonds(ref))) {
            Bond *nbrbond = ref->getOwningMol()[nbrbondItr];
            if (&*(nbrbond) == &*bond || !nbrbond->getIsAromatic())
              continue;
            Atom *nbrnbr = nbrbond->getOtherAtom(ref);
            unsigned int nbrnbrIdx = nbrnbr->getIdx();
            if (!needs_dbl_bond->BitIsSet(nbrnbrIdx))
              continue;
            degrees[nbrnbrIdx]--;
            if (degrees[nbrnbrIdx] == 1) {
              degreeOneAtoms.push_back(nbrnbr);
              change = true;
            }
          }
        }
        // only a single double bond can be made to atom so we can break here
        break;
      }
      if (change)
        break;
    }

    // We exit if we are finished or if no degree 2/3 nodes can be set
    if (!change)
      break;
  }

  // Tidy up
  free(degrees);

  return finished;
}

// The isDoubleBond alternates between double and single, as we need to find
// an alternating path
bool Kekulizer::FindPath(unsigned int atomidx, bool isDoubleBond,
                         OBBitVec &visited) {
  if (needs_dbl_bond->BitIsSet(atomidx))
    return true;
  visited.SetBitOn(atomidx);
  Atom *atom = m_mol->getAtomWithIdx(atomidx);
  for (const auto &bondItr :
       boost::make_iterator_range(atom->getOwningMol().getAtomBonds(atom))) {
    Bond *bond = atom->getOwningMol()[bondItr];
    if (!bond->getIsAromatic())
      continue;
    Atom *nbr = bond->getOtherAtom(atom);
    if (!kekule_system->BitIsSet(nbr->getIdx()))
      continue;
    if (doubleBonds->BitIsSet(bond->getIdx()) == isDoubleBond) {
      if (visited.BitIsSet(nbr->getIdx()))
        continue;
      bool found_path = FindPath(nbr->getIdx(), !isDoubleBond, visited);
      if (found_path) {
        m_path.push_back(nbr->getIdx());
        return true;
      }
    }
  }
  visited.SetBitOff(atomidx);
  return false;
}

bool Kekulizer::BackTrack() {
  // With an odd number of bits, it's never going to kekulize fully, but let's
  // fill in as many as we can
  unsigned int count = needs_dbl_bond->CountBits();

  unsigned int total_handled = 0;
  int idx;
  for (idx = needs_dbl_bond->FirstBit(); idx != needs_dbl_bond->EndBit();
       idx = needs_dbl_bond->NextBit(idx)) {
    total_handled++;
    // If there is no additional bit available to match this bit, then terminate
    if (total_handled == count)
      return false;

    // Our goal is to find an alternating path to another atom
    // that needs a double bond
    needs_dbl_bond->SetBitOff(
        idx); // to avoid the trivial null path being found
    OBBitVec visited(atomArraySize);
    m_path.clear();
    bool found_path = FindPath(idx, false, visited);
    if (!found_path) {               // could only happen if not kekulizable
      needs_dbl_bond->SetBitOn(idx); // reset
      continue;
    }
    total_handled++;
    m_path.push_back(idx);
    needs_dbl_bond->SetBitOff(m_path[0]);
    // Flip all of the bond orders on the path from double<-->single
    for (unsigned int i = 0; i < m_path.size() - 1; ++i) {
      Bond *bond = m_mol->getBondBetweenAtoms(m_path[i], m_path[i + 1]);
      if (i % 2 == 0)
        doubleBonds->SetBitOn(bond->getIdx());
      else
        doubleBonds->SetBitOff(bond->getIdx());
    }
  }
  return needs_dbl_bond->IsEmpty();
}

// I'd like to thank John Mayfield for many helpful discussions on the topic of
// kekulization, without which this implementation would not have been
// possible.
//
// OBKekulize() implements a two-step kekulization:
//   Step one: try a greedy match
//   Step two: try an exhaustive backtracking (using the results of step one)
//
// The greedy match algorithm is outlined in the thesis of John May
// and indeed NeedsDoubleBond() is based on the implementation in Beam.
// The greedy algorithm almost always works. But when it doesn't, step two is
// needed.
//
// The goal of the exhaustive backtracking is find a path of alternating
// single/double bonds between two radicals and flip those bonds. For more
// information, read about augmenting paths in the context of perfect matching.
// John's thesis instead describes the use of Edmond's Blossom algorithm which
// scales better - this may or may not be faster in practice for typical
// chemical graphs.
//
// Potential speedups:
//   * Is OBBitVec performant? I don't know - it seems to do a lot of bounds
//   checking.
//     You could try replacing all usages theoreof with
//     std::vector<char>, where the char could possibly handle several flags.
//   * Before trying the exhaustive search, try a BFS. I have a feeling that
//   this would work
//     90% of the time.
//   * There's a lot of switching between atoms and atom indices (and similar
//   for bonds).
//     Was this completely necessary?
//   * The iterator over degree 2 and 3 nodes may iterate twice - it would have
//   been
//     faster if I just took the first degree 2 or 3 node I came across, but
//     would it have worked as well?

bool OBKekulize(RWMol *mol) {
  Kekulizer kekulizer(mol);

  bool success = kekulizer.GreedyMatch();
  if (!success) {
    success = kekulizer.BackTrack();
  }

  kekulizer.AssignDoubleBonds();
  return success;
}

//! \file kekulize.cpp
//! \brief Algorithm to kekulize a molecule
