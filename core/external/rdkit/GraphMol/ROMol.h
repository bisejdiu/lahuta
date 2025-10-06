//
//  Copyright (C) 2003-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file ROMol.h

  \brief Defines the primary molecule class \c ROMol as well as associated
  typedefs

*/

// ----------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------
// This file has been heavily modified from the original to fit with Lahuta's, highly speicific objectives.
// The original implementation relied on an adjacency list graph to store atom and bond information.
// The flexibility that this afforded was not needed, and the lack of control over how the graph was
// constructed was frustrating. I have, in addition to the adjacency list, also added a CSR graph.
// This makes graph construction easier and faster. The redesign has been anything but easy, and the
// resulting code is not as clean as I would have liked. To keep changes contained, I had to avoid
// templating ROMol. Iterators, in particular, were a mess to deal with.
//
// A few caveats are in order:
// - The code has been adapted to fit with our intended use of the library. While I have attempted to maintain
//   the original functionality, I have not tested all of the features, and I would not be surprised if some
//   of them are broken. They will be fixed if and when they are needed.
// - Old style iterations (e.g. beginAtoms(), endAtoms()) is likely broken and should be avoided. In any case,
//   there really is not much need for them anymore. Iterators that allowed any type of iteration over atoms:
//   e.g.: for (auto &atom : mol.atoms()) { ... }; or for (auto *atom : mol.atoms()) { ... }; won't work
//   anymore. You get a pointer to the object and you should be happy with that.
// - getAtomBonds() and getAtomNeighbors() are a bit weird in that they return a pair of hardcoded iterators.
//   I tried supporting them but that lead to massive complications either downstream or upstream. I regret
//   even trying. In any case, these are not really striclty needed, given that they serve the same purpose as
//   atomBonds and atomNeighbors. For now they work with the default graph type, but at some point they should
//   be removed.
// - Finally, it would be nice to have conversion functions between the two graph types. But given the special
//   purpose redesign, I don't see the point. If the need arises, we can implement them later. For now, it
//   would be counterproductive to try to add features for a special purpose redesign to fit a specific use
//   case.
// - While we have direct control over graph construction, this comes at the cost of pointer indirection (via
//   ptr to impl) and vtable lookups. In practice, these are quite cheap. The cost is more on mental overhead.
//                                                                                      -- Besian, April 2025
// ----------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------

#ifndef RD_ROMOL_H
#define RD_ROMOL_H

/// Std stuff
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <cstddef>
#include <utility>
#include <map>

// boost stuff
#include <RDGeneral/BoostStartInclude.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/dynamic_bitset.hpp>

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <boost/serialization/split_member.hpp>
#endif
#include <RDGeneral/BoostEndInclude.h>

// our stuff
#include <RDGeneral/types.h>
#include <RDGeneral/RDProps.h>
#include "Atom.h"
#include "Bond.h"
#include "Conformer.h"
#include "SubstanceGroup.h"
#include "StereoGroup.h"
#include "RingInfo.h"
#include "Graphs.hpp"

namespace RDKit {
class SubstanceGroup;
class Atom;
class Bond;
class MolPickler;
class RWMol;
class QueryAtom;
class QueryBond;
class RingInfo;

template <class T1, class T2>
class AtomIterator_;
class BondIterator_;
class ConstBondIterator_;

template <class T1, class T2>
class AromaticAtomIterator_;
template <class T1, class T2>
class HeteroatomIterator_;
template <class T1, class T2>
class QueryAtomIterator_;
template <class T1, class T2>
class MatchingAtomIterator_;

const int ci_RIGHTMOST_ATOM = -0xBADBEEF;
const int ci_LEADING_BOND   = -0xBADBEEF + 1;
const int ci_ATOM_HOLDER    = -0xDEADD06;

class ROMol : public RDProps {
public:
  friend class MolPickler;
  friend class RWMol;

  using vertex_descriptor = typename boost::graph_traits<MolGraph>::vertex_descriptor;
  using edge_descriptor   = typename boost::graph_traits<MolGraph>::edge_descriptor;

  using ADJ_ITER    = typename boost::graph_traits<MolGraph>::adjacency_iterator;
  using EDGE_ITER   = typename boost::graph_traits<MolGraph>::edge_iterator;
  using OEDGE_ITER  = typename boost::graph_traits<MolGraph>::out_edge_iterator;
  using VERTEX_ITER = typename boost::graph_traits<MolGraph>::vertex_iterator;

  using BOND_ITER_PAIR  = std::pair<EDGE_ITER, EDGE_ITER>;
  using ATOM_ITER_PAIR  = std::pair<VERTEX_ITER, VERTEX_ITER>;
  using OBOND_ITER_PAIR = std::pair<OEDGE_ITER, OEDGE_ITER>;
  using ADJ_ITER_PAIR   = std::pair<ADJ_ITER, ADJ_ITER>;

  using ATOM_PTR_VECT = std::vector<Atom *>;
  using BOND_PTR_VECT = std::vector<Bond *>;
  using ATOM_PTR_LIST = std::list<Atom *>;
  using BOND_PTR_LIST = std::list<Bond *>;
  using ATOM_BOOKMARK_MAP = std::map<int, ATOM_PTR_LIST>;
  using BOND_BOOKMARK_MAP = std::map<int, BOND_PTR_LIST>;
  using CONF_SPTR_LIST = std::list<CONFORMER_SPTR>;

  typedef class AtomIterator_<Atom, ROMol> AtomIterator;
  typedef class AtomIterator_<const Atom, const ROMol> ConstAtomIterator;
  typedef class BondIterator_ BondIterator;
  typedef class ConstBondIterator_ ConstBondIterator;

  typedef CONF_SPTR_LIST::iterator CONF_SPTR_LIST_I;
  typedef CONF_SPTR_LIST::const_iterator CONF_SPTR_LIST_CI;
  typedef std::pair<CONF_SPTR_LIST_I, CONF_SPTR_LIST_I> CONFS_I_PAIR;

  typedef CONF_SPTR_LIST_I ConformerIterator;
  typedef CONF_SPTR_LIST_CI ConstConformerIterator;

  //! get an AtomIterator pointing at our first Atom
  AtomIterator beginAtoms();
  ConstAtomIterator beginAtoms() const;
  //! get an AtomIterator pointing at the end of our Atoms
  AtomIterator endAtoms();
  ConstAtomIterator endAtoms() const;

  //! get a BondIterator pointing at our first Bond
  BondIterator beginBonds();
  ConstBondIterator beginBonds() const;
  //! get a BondIterator pointing at the end of our Bonds
  BondIterator endBonds();
  ConstBondIterator endBonds() const;

  //! get an ConformerIterator pointing at our first Conformer
  ConformerIterator beginConformers() { return d_confs.begin(); }
  ConformerIterator endConformers() { return d_confs.end(); }

  //! get a ConstConformerIterator pointing at our first Conformer
  ConstConformerIterator beginConformers() const { return d_confs.begin(); }
  ConstConformerIterator endConformers() const { return d_confs.end(); }

  void skipAtomCleanupDespiteOwnership(bool skip) {
    should_delete_atoms = !skip;
  }
  void skipBondCleanupDespiteOwnership(bool skip) {
    should_delete_bonds = !skip;
  }

  ROMol() : RDProps(), m_impl(std::make_unique<MolGraphImpl>()) { initMol(); }

  explicit ROMol(GraphType type) {
    switch (type) {
      case GraphType::MolGraph:
        m_impl = std::make_unique<MolGraphImpl>();
        break;
      case GraphType::CSRMolGraph:
        m_impl = std::make_unique<CSRMolGraphImpl>();
        break;
    }
    initMol();
  }

  // NOTE: Will not take ownership of the atoms or bonds
  ROMol(
      const std::vector<Atom *> &atoms,
      const std::vector<Bond *> &bonds,
      GraphType type = GraphType::CSRMolGraph) {

    initMol();

    for (auto atom : atoms) atom->setOwningMol(this);
    for (auto bond : bonds) bond->setOwningMol(this);

    switch (type) {
      case GraphType::MolGraph: {
        std::vector<std::tuple<size_t, size_t, Bond *>> bonds_list;
        bonds_list.reserve(bonds.size());
        for (auto bond : bonds) {
          bonds_list.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx(), bond);
        }

        m_impl = std::make_unique<MolGraphImpl>();
        auto *impl = static_cast<MolGraphImpl *>(m_impl.get());
        impl->build(atoms, bonds_list);
        break;
      }
      case GraphType::CSRMolGraph: {
        m_impl = std::make_unique<CSRMolGraphImpl>();
        auto *impl = static_cast<CSRMolGraphImpl *>(m_impl.get());
        impl->build(atoms, bonds);
        break;
      }
    }
    should_delete_bonds = false;
    should_delete_atoms = false;
  }

  ROMol(const ROMol &other, bool quickCopy = false, int confId = -1)
      : RDProps(), m_impl(std::make_unique<MolGraphImpl>()) {
    initFromOther(other, quickCopy, confId);
    numBonds = m_impl->getNumEdges();
  }

  ROMol(const std::string &binStr) : RDProps(), m_impl(std::make_unique<MolGraphImpl>()) {
    initMol();
    numBonds = m_impl->getNumEdges();
  }

  ROMol(const std::string &binStr, unsigned int propertyFlags)
      : RDProps(), m_impl(std::make_unique<MolGraphImpl>()) {
    initMol();
    numBonds = m_impl->getNumEdges();
  }

  ROMol(ROMol &&o) noexcept : RDProps(std::move(o)) {
    m_impl = std::move(o.m_impl);
    should_delete_atoms = o.should_delete_atoms;
    should_delete_bonds = o.should_delete_bonds;
    d_atomBookmarks = std::move(o.d_atomBookmarks);
    d_bondBookmarks = std::move(o.d_bondBookmarks);
    d_confs = std::move(o.d_confs);
    d_sgroups = std::move(o.d_sgroups);
    d_stereo_groups = std::move(o.d_stereo_groups);
    numBonds = o.numBonds;
    o.numBonds = 0;
    dp_ringInfo = std::exchange(o.dp_ringInfo, nullptr);
    dp_delAtoms = std::exchange(o.dp_delAtoms, nullptr);
    dp_delBonds = std::exchange(o.dp_delBonds, nullptr);

    for (auto atom : atoms()) {
      atom->setOwningMol(this);
    }

    for (auto bond : bonds()) {
      bond->setOwningMol(this);
    }

    for (auto conf : d_confs) {
      conf->setOwningMol(this);
    }

    for (auto &sg : d_sgroups) {
      sg.setOwningMol(this);
    }
  }

  ROMol &operator=(ROMol &&o) noexcept {
    if (this == &o) return *this;

    RDProps::operator=(std::move(o));
    m_impl = std::move(o.m_impl);
    should_delete_atoms = o.should_delete_atoms;
    should_delete_bonds = o.should_delete_bonds;
    d_atomBookmarks = std::move(o.d_atomBookmarks);
    d_bondBookmarks = std::move(o.d_bondBookmarks);
    d_confs = std::move(o.d_confs);
    d_sgroups = std::move(o.d_sgroups);
    d_stereo_groups = std::move(o.d_stereo_groups);
    numBonds = o.numBonds;
    o.numBonds = 0;
    dp_ringInfo = std::exchange(o.dp_ringInfo, nullptr);
    dp_delAtoms = std::exchange(o.dp_delAtoms, nullptr);
    dp_delBonds = std::exchange(o.dp_delBonds, nullptr);

    for (auto atom : atoms()) {
      atom->setOwningMol(this);
    }

    for (auto bond : bonds()) {
      bond->setOwningMol(this);
    }

    for (auto conf : d_confs) {
      conf->setOwningMol(this);
    }

    for (auto &sg : d_sgroups) {
      sg.setOwningMol(this);
    }

    return *this;
  }

  ROMol &operator=(const ROMol &) = delete;
  virtual ~ROMol() { ROMol::destroy(); }

  unsigned int getNumAtoms() const { return m_impl->getNumVertices(); }

  // FIX: I'm not sure we handle this correctly!
  unsigned int getNumAtoms(bool onlyExplicit) const {
    auto res = m_impl->getNumVertices();
    if (!onlyExplicit) {
      for (const auto atom : atoms()) {
        res += atom->getTotalNumHs();
      }
    }
    return res;
  }

  unsigned int getNumHeavyAtoms() const {
    unsigned int res = 0;
    for (const auto atom : atoms()) {
      if (atom->getAtomicNum() > 1) ++res;
    }
    return res;
  }

  Atom *getAtomWithIdx(unsigned int idx) {
    URANGE_CHECK(idx, getNumAtoms());
    auto res = m_impl->getAtom(idx);
    POSTCONDITION(res, "");
    return res;
  }

  const Atom *getAtomWithIdx(unsigned int idx) const {
    URANGE_CHECK(idx, getNumAtoms());
    auto res = m_impl->getAtom(idx);
    POSTCONDITION(res, "");
    return res;
  }

  unsigned int getAtomDegree(const Atom *at) const {
    PRECONDITION(at, "no atom");
    PRECONDITION(&at->getOwningMol() == this, "atom not associated with this molecule");
    return m_impl->getAtomDegree(at->getIdx());
  }

  unsigned int getNumBonds(bool onlyHeavy = true) const {
    // By default return the bonds that connect only the heavy atoms
    // hydrogen connecting bonds are ignores
    unsigned int res = m_impl->getNumEdges();
    if (!onlyHeavy) {
      for (auto atom : atoms()) {
        res += atom->getTotalNumHs();
      }
    }
    return res;
  }

  Bond *getBondWithIdx(unsigned int idx) {
    URANGE_CHECK(idx, getNumBonds());
    auto res = m_impl->getBond(idx);
    POSTCONDITION(res != nullptr, "Invalid bond requested");
    return res;
  }

  const Bond *getBondWithIdx(unsigned int idx) const {
    URANGE_CHECK(idx, getNumBonds());
    auto res = m_impl->getBond(idx);
    POSTCONDITION(res != nullptr, "Invalid bond requested");
    return res;
  }

  Bond *getBondBetweenAtoms(unsigned int idx1, unsigned int idx2) {
    URANGE_CHECK(idx1, getNumAtoms());
    URANGE_CHECK(idx2, getNumAtoms());
    auto [e, success] = m_impl->getBondBetweenAtoms(idx1, idx2);
    return success ? static_cast<Bond *>(e) : nullptr;
  }

  const Bond *getBondBetweenAtoms(unsigned int idx1, unsigned int idx2) const {
    URANGE_CHECK(idx1, getNumAtoms());
    URANGE_CHECK(idx2, getNumAtoms());
    auto [e, success] = m_impl->getBondBetweenAtoms(idx1, idx2);
    return success ? static_cast<const Bond *>(e) : nullptr;
  }

  int getConnectedComponents(std::vector<int> &mapping) const {
    return m_impl->getConnectedComponents(mapping);
  }

  void setAtomBookmark(Atom *at, int mark) { d_atomBookmarks[mark].push_back(at); }
  void replaceAtomBookmark(Atom *at, int mark) {
    d_atomBookmarks[mark].clear();
    d_atomBookmarks[mark].push_back(at);
  }

  Atom *getAtomWithBookmark(int mark) {
    auto it = d_atomBookmarks.find(mark);
    PRECONDITION((it != d_atomBookmarks.end() && !it->second.empty()), "atom bookmark not found");
    return it->second.front();
  }

  Atom *getUniqueAtomWithBookmark(int mark) {
    auto it = d_atomBookmarks.find(mark);
    PRECONDITION(it != d_atomBookmarks.end(), "bookmark not found");
    PRECONDITION(it->second.size() == 1, "unique atom bookmark not found");
    return it->second.front();
  }

  ATOM_PTR_LIST &getAllAtomsWithBookmark(int mark) {
    auto it = d_atomBookmarks.find(mark);
    PRECONDITION(it != d_atomBookmarks.end(), "atom bookmark not found");
    return it->second;
  }

  void clearAtomBookmark(int mark) { d_atomBookmarks.erase(mark); }

  void clearAtomBookmark(int mark, const Atom *atom) {
    PRECONDITION(atom, "no atom");
    auto it = d_atomBookmarks.find(mark);
    if (it != d_atomBookmarks.end()) {
      auto &lst = it->second;
      unsigned int tgtIdx = atom->getIdx();
      auto entry =
          std::find_if(lst.begin(), lst.end(), [&tgtIdx](auto ptr) { return ptr->getIdx() == tgtIdx; });
      if (entry != lst.end()) {
        lst.erase(entry);
      }
      if (lst.empty()) {
        d_atomBookmarks.erase(mark);
      }
    }
  }
  void clearAllAtomBookmarks() { d_atomBookmarks.clear(); }
  bool hasAtomBookmark(int mark) const { return (d_atomBookmarks.count(mark) != 0); }
  ATOM_BOOKMARK_MAP *getAtomBookmarks() { return &d_atomBookmarks; }

  void setBondBookmark(Bond *bond, int mark) { d_bondBookmarks[mark].push_back(bond); }

  Bond *getBondWithBookmark(int mark) {
    auto it = d_bondBookmarks.find(mark);
    PRECONDITION(it != d_bondBookmarks.end() && !it->second.empty(), "bond bookmark not found");
    return it->second.front();
  }

  Bond *getUniqueBondWithBookmark(int mark) {
    auto it = d_bondBookmarks.find(mark);
    PRECONDITION(it != d_bondBookmarks.end(), "bookmark not found");
    PRECONDITION(it->second.size() == 1, "unique bond bookmark not found");
    return it->second.front();
  }

  BOND_PTR_LIST &getAllBondsWithBookmark(int mark) {
    auto it = d_bondBookmarks.find(mark);
    PRECONDITION(it != d_bondBookmarks.end(), "bond bookmark not found");
    return it->second;
  }

  void clearBondBookmark(int mark) { d_bondBookmarks.erase(mark); }

  void clearBondBookmark(int mark, const Bond *bond) {
    PRECONDITION(bond, "no bond");
    auto it = d_bondBookmarks.find(mark);
    if (it != d_bondBookmarks.end()) {
      auto &lst = it->second;
      unsigned int tgtIdx = bond->getIdx();
      auto entry =
          std::find_if(lst.begin(), lst.end(), [&tgtIdx](auto ptr) { return ptr->getIdx() == tgtIdx; });
      if (entry != lst.end()) {
        lst.erase(entry);
      }
      if (lst.empty()) {
        d_bondBookmarks.erase(mark);
      }
    }
  }

  void clearAllBondBookmarks() { d_bondBookmarks.clear(); }
  bool hasBondBookmark(int mark) const { return (d_bondBookmarks.count(mark) != 0); }
  BOND_BOOKMARK_MAP *getBondBookmarks() { return &d_bondBookmarks; }

  const Conformer &getConformer(int id = -1) const {
    // make sure we have more than one conformation
    if (d_confs.size() == 0) {
      throw ConformerException("No conformations available on the molecule");
    }

    if (id < 0) {
      return *(d_confs.front());
    }
    auto cid = (unsigned int)id;
    for (auto conf : d_confs) {
      if (conf->getId() == cid) {
        return *conf;
      }
    }
    // we did not find a conformation with the specified ID
    std::string mesg = "Can't find conformation with ID: ";
    mesg += id;
    throw ConformerException(mesg);
  }
  Conformer &getConformer(int id = -1) {
    if (d_confs.size() == 0) {
      throw ConformerException("No conformations available on the molecule");
    }

    if (id < 0) {
      return *(d_confs.front());
    }
    auto cid = (unsigned int)id;
    for (auto conf : d_confs) {
      if (conf->getId() == cid) {
        return *conf;
      }
    }
    // we did not find a conformation with the specified ID
    std::string mesg = "Can't find conformation with ID: ";
    mesg += id;
    throw ConformerException(mesg);
  }
  void removeConformer(unsigned int id) {
    for (auto ci = d_confs.begin(); ci != d_confs.end(); ++ci) {
      if ((*ci)->getId() == id) {
        d_confs.erase(ci);
        return;
      }
    }
  }
  void clearConformers() { d_confs.clear(); }
  unsigned int addConformer(Conformer *conf, bool assignId = false) {
    PRECONDITION(conf, "bad conformer");
    PRECONDITION(conf->getNumAtoms() == this->getNumAtoms(), "Number of atom mismatch");
    if (assignId) {
      int maxId = -1;
      for (auto cptr : d_confs) {
        maxId = std::max((int)(cptr->getId()), maxId);
      }
      maxId++;
      conf->setId((unsigned int)maxId);
    }
    conf->setOwningMol(static_cast<ROMol *>(this));
    CONFORMER_SPTR nConf(conf);
    d_confs.push_back(nConf);
    return conf->getId();
  }
  unsigned int getNumConformers() const { return (unsigned int)d_confs.size(); }

  RingInfo *getRingInfo() const { return dp_ringInfo; }

  std::pair<ADJ_ITER, ADJ_ITER> getAtomNeighbors(const Atom *at) const {
    if (auto *molGraphImpl = static_cast<MolGraphImpl *>(m_impl.get())) {
      return molGraphImpl->getMolGraphAtomNeighbors(at->getIdx());
    } else if (auto *csrGraphImpl = static_cast<CSRMolGraphImpl *>(m_impl.get())) {
      // we need to throw because we can't (don't want to) convert iterators
      throw std::runtime_error("CSR graph not supported in getAtomNeighbors");
    }
    throw std::runtime_error("Unknown graph implementation");
  }

  std::pair<OEDGE_ITER, OEDGE_ITER> getAtomBonds(const Atom *at) const {
    unsigned int idx = static_cast<unsigned int>(at->getIdx());

    if (auto *molGraphImpl = static_cast<MolGraphImpl *>(m_impl.get())) {
      return molGraphImpl->getMolGraphAtomBonds(idx);
    } else if (auto *csrGraphImpl = static_cast<CSRMolGraphImpl *>(m_impl.get())) {
      // we need to throw because we can't (don't want to) convert iterators
      throw std::runtime_error("CSR graph not supported in getAtomBonds");
    } else {
      throw std::runtime_error("Unknown graph implementation");
    }
  }

  BOND_ITER_PAIR getEdges() {
    if (auto *molGraphImpl = static_cast<MolGraphImpl *>(m_impl.get())) {
      return boost::edges(molGraphImpl->getGraph());
    }
    throw std::runtime_error("Unknown graph implementation");
  }

  BOND_ITER_PAIR getEdges() const {
    if (auto *molGraphImpl = static_cast<MolGraphImpl *>(m_impl.get())) {
      return boost::edges(molGraphImpl->getGraph());
    }
    throw std::runtime_error("Unknown graph implementation");
  }

  const MolGraph &getTopology() const {
    if (auto *molGraphImpl = static_cast<MolGraphImpl *>(m_impl.get())) {
      return molGraphImpl->getGraph();
    }
    throw std::runtime_error("Unknown graph implementation");
  }

  bool hasQuery() const {
    for (auto atom : atoms()) {
      if (atom->hasQuery()) {
        return true;
      }
    }
    for (auto bond : bonds()) {
      if (bond->hasQuery()) {
        return true;
      }
    }
    return false;
  }

  lt::AtomRange atoms() const {
    auto [begin, end] = m_impl->getAtomIterators();
    return {std::move(begin), std::move(end)};
  }

  lt::BondRange bonds() const {
    auto [begin, end] = m_impl->getBondIterators();
    return {std::move(begin), std::move(end)};
  }

  lt::BondRange atomBonds(const Atom *atom) const {
    auto [begin, end] = m_impl->getAtomBonds(atom->getIdx());
    return {std::move(begin), std::move(end)};
  }

  lt::AtomRange atomNeighbors(const Atom *atom) const {
    auto [begin, end] = m_impl->getAtomNeighbors(atom->getIdx());
    return {std::move(begin), std::move(end)};
  }

  void clearComputedProps(bool includeRings = true) const {
    // the SSSR information:
    if (includeRings) {
      this->dp_ringInfo->reset();
    }

    RDProps::clearComputedProps();

    for (auto atom : atoms()) {
      atom->getProps()->clearComputedProps();
    }

    for (auto bond : bonds()) {
      bond->getProps()->clearComputedProps();
    }
  }
  void updatePropertyCache(bool strict = true) {
    for (auto atom : atoms()) {
      atom->updatePropertyCache(strict);
    }
    for (auto bond : bonds()) {
      bond->updatePropertyCache(strict);
    }
  }
  bool needsUpdatePropertyCache() const {
    for (const auto atom : atoms()) {
      if (atom->needsUpdatePropertyCache()) {
        return true;
      }
    }
    // there is no test for bonds yet since they do not obtain a valence property
    return false;
  }

  void debugMol(std::ostream &str) const {
    str << "Atoms:" << std::endl;
    for (const auto atom : atoms()) {
      str << "\t" << *atom << std::endl;
    }

    str << "Bonds:" << std::endl;
    for (const auto bond : bonds()) {
      str << "\t" << *bond << std::endl;
    }

    const auto &sgs = getSubstanceGroups();
    if (!sgs.empty()) {
      str << "Substance Groups:" << std::endl;
      for (const auto &sg : sgs) {
        str << "\t" << sg << std::endl;
      }
    }

    const auto &stgs = getStereoGroups();
    if (!stgs.empty()) {
      unsigned idx = 0;
      str << "Stereo Groups:" << std::endl;
      for (const auto &stg : stgs) {
        str << "\t" << idx << ' ' << stg << std::endl;
        ++idx;
      }
    }
  }

  Atom *operator[](const vertex_descriptor &v) { return m_impl->getAtom(v); }
  const Atom *operator[](const vertex_descriptor &v) const { return m_impl->getAtom(v); }

  Bond *operator[](const edge_descriptor &e) {
    if (m_impl->isMolGraph()) {
      return m_impl->getBond(e);
    }
    auto src = boost::source(e, MolGraph());
    auto dst = boost::target(e, MolGraph());

    auto [bond, exists] = m_impl->getBondBetweenAtoms(src, dst);
    if (!exists) {
      throw std::runtime_error("Bond not found");
    }
    return bond;
  }

  const Bond *operator[](const edge_descriptor &e) const {
    if (m_impl->isMolGraph()) {
      return m_impl->getBond(e);
    }
    auto src = boost::source(e, MolGraph());
    auto dst = boost::target(e, MolGraph());

    auto [bond, exists] = m_impl->getBondBetweenAtoms(src, dst);
    if (!exists) {
      throw std::runtime_error("Bond not found");
    }
    return bond;
  }

  const std::vector<StereoGroup> &getStereoGroups() const { return d_stereo_groups; }
  void setStereoGroups(std::vector<StereoGroup> sgs) { d_stereo_groups = std::move(sgs); }

  const std::vector<SubstanceGroup> &getSubstanceGroups() const { return d_sgroups; }
  std::vector<SubstanceGroup> &getSubstanceGroups() { return d_sgroups; }

  GraphImpl *getImpl() { return m_impl.get(); }
  const GraphImpl *getImpl() const { return m_impl.get(); }

private:
  std::unique_ptr<GraphImpl> m_impl;
  bool should_delete_atoms{true};
  bool should_delete_bonds{true};

  ATOM_BOOKMARK_MAP d_atomBookmarks;
  BOND_BOOKMARK_MAP d_bondBookmarks;
  RingInfo *dp_ringInfo = nullptr;
  CONF_SPTR_LIST d_confs;
  std::vector<SubstanceGroup> d_sgroups;
  std::vector<StereoGroup> d_stereo_groups;
  std::unique_ptr<boost::dynamic_bitset<>> dp_delAtoms = nullptr;
  std::unique_ptr<boost::dynamic_bitset<>> dp_delBonds = nullptr;
  friend std::vector<SubstanceGroup> &getSubstanceGroups(ROMol &);
  friend const std::vector<SubstanceGroup> &getSubstanceGroups(const ROMol &);

protected:
  unsigned int numBonds{0};

  void initMol() {
    d_props.reset();
    dp_ringInfo = new RingInfo();
    // ok every molecule contains a property entry called
    // RDKit::detail::computedPropName
    // which provides
    //  list of property keys that correspond to value that have been computed
    // this can used to blow out all computed properties while leaving the rest
    // along
    // initialize this list to an empty vector of strings
    STR_VECT computed;
    d_props.setVal(RDKit::detail::computedPropName, computed);
  }

  virtual void destroy() {

    d_atomBookmarks.clear();
    d_bondBookmarks.clear();

    if (m_impl && should_delete_bonds) {
      for (auto bond : bonds()) {
        delete bond;
      }
    }

    if (m_impl && should_delete_atoms) {
      for (auto atom : atoms()) {
        delete atom;
      }
    }

    if (dp_ringInfo) {
      delete dp_ringInfo;
      dp_ringInfo = nullptr;
    }

    d_sgroups.clear();
    d_stereo_groups.clear();

    m_impl = nullptr;
  }

  void initFromOther(const ROMol &other, bool quickCopy, int confId) {
    if (this == &other) {
      return;
    }
    numBonds = 0;

    // copy over the atoms
    for (const auto oatom : other.atoms()) {
      constexpr bool updateLabel = false;
      constexpr bool takeOwnership = true;
      addAtom(oatom->copy(), updateLabel, takeOwnership);
    }

    // and the bonds:
    for (const auto obond : other.bonds()) {
      addBond(obond->copy(), true);
    }

    // ring information
    delete dp_ringInfo;
    if (other.dp_ringInfo) {
      dp_ringInfo = new RingInfo(*(other.dp_ringInfo));
    } else {
      dp_ringInfo = new RingInfo();
    }

    // enhanced stereochemical information
    d_stereo_groups.clear();
    for (auto &otherGroup : other.d_stereo_groups) {
      std::vector<Atom *> atoms;
      for (auto &otherAtom : otherGroup.getAtoms()) {
        atoms.push_back(getAtomWithIdx(otherAtom->getIdx()));
      }
      std::vector<Bond *> bonds;
      for (auto &otherBond : otherGroup.getBonds()) {
        bonds.push_back(getBondWithIdx(otherBond->getIdx()));
      }
      d_stereo_groups.emplace_back(
          otherGroup.getGroupType(),
          std::move(atoms),
          std::move(bonds),
          otherGroup.getReadId());
      d_stereo_groups.back().setWriteId(otherGroup.getWriteId());
    }

    if (other.dp_delAtoms) {
      dp_delAtoms.reset(new boost::dynamic_bitset<>(*other.dp_delAtoms));
    } else {
      dp_delAtoms.reset(nullptr);
    }
    if (other.dp_delBonds) {
      dp_delBonds.reset(new boost::dynamic_bitset<>(*other.dp_delBonds));
    } else {
      dp_delBonds.reset(nullptr);
    }

    if (!quickCopy) {
      // copy conformations
      for (const auto &conf : other.d_confs) {
        if (confId < 0 || rdcast<int>(conf->getId()) == confId) {
          this->addConformer(new Conformer(*conf));
        }
      }

      // Copy sgroups
      for (const auto &sg : getSubstanceGroups()) {
        addSubstanceGroup(*this, sg);
      }

      d_props = other.d_props;

      // Bookmarks should be copied as well:
      for (auto abmI : other.d_atomBookmarks) {
        for (const auto *aptr : abmI.second) {
          setAtomBookmark(getAtomWithIdx(aptr->getIdx()), abmI.first);
        }
      }
      for (auto bbmI : other.d_bondBookmarks) {
        for (const auto *bptr : bbmI.second) {
          setBondBookmark(getBondWithIdx(bptr->getIdx()), bbmI.first);
        }
      }
    } else {
      d_props.reset();
      STR_VECT computed;
      d_props.setVal(RDKit::detail::computedPropName, computed);
    }
    // std::cerr<<"---------    done init from other: "<<this<<"
    // "<<&other<<std::endl;
  }

  unsigned int addAtom(Atom *atom_pin, bool updateLabel, bool takeOwnership = false) {
    PRECONDITION(atom_pin, "null atom passed in");
    PRECONDITION(!takeOwnership || !atom_pin->hasOwningMol() ||
                     &atom_pin->getOwningMol() == this,
                 "cannot take ownership of an atom which already has an owner");
    Atom *atom_p;
    if (!takeOwnership) {
      atom_p = atom_pin->copy();
    } else {
      atom_p = atom_pin;
    }

    atom_p->setOwningMol(this);

    if (!m_impl) throw ValueErrorException("Graph implementation not set");
    auto which = m_impl->addAtom(atom_p);

    atom_p->setIdx(which);
    if (updateLabel) {
      replaceAtomBookmark(atom_p, ci_RIGHTMOST_ATOM);
    }
    for (auto &conf : d_confs) {
      conf->setAtomPos(which, RDGeom::Point3D(0.0, 0.0, 0.0));
    }
    return rdcast<unsigned int>(which);
  };

  unsigned int addBond(Bond *bond_pin, bool takeOwnership = false) {
    PRECONDITION(bond_pin, "null bond passed in");
    PRECONDITION(!takeOwnership || !bond_pin->hasOwningMol() ||
                     &bond_pin->getOwningMol() == this,
                 "cannot take ownership of an bond which already has an owner");
    URANGE_CHECK(bond_pin->getBeginAtomIdx(), getNumAtoms());
    URANGE_CHECK(bond_pin->getEndAtomIdx(), getNumAtoms());
    PRECONDITION(bond_pin->getBeginAtomIdx() != bond_pin->getEndAtomIdx(), "attempt to add self-bond");
    PRECONDITION(!m_impl->getBondBetweenAtoms(bond_pin->getBeginAtomIdx(), bond_pin->getEndAtomIdx()).second,
                 "bond already exists");

    Bond *bond_p;
    if (!takeOwnership) {
      bond_p = bond_pin->copy();
    } else {
      bond_p = bond_pin;
    }
    bond_p->setOwningMol(this);
    m_impl->addBond(bond_p->getBeginAtomIdx(), bond_p->getEndAtomIdx(), bond_p);
    bond_p->setIdx(numBonds);
    numBonds++;
    return numBonds;
  }

  // all checks are removed. you should know what you are doing
  unsigned int addBonds(std::vector<Bond *> &bonds) {
    for (auto bond_p : bonds) {
      bond_p->setOwningMol(this);
      m_impl->addBond(bond_p->getBeginAtomIdx(), bond_p->getEndAtomIdx(), bond_p);
      bond_p->setIdx(numBonds);
      numBonds++;
    }

    should_delete_bonds = false;
    return numBonds;
  }
};

typedef std::vector<ROMol> MOL_VECT;
typedef boost::shared_ptr<ROMol> ROMOL_SPTR;
typedef std::vector<ROMol *> MOL_PTR_VECT;
typedef std::vector<ROMOL_SPTR> MOL_SPTR_VECT;

typedef MOL_PTR_VECT::const_iterator MOL_PTR_VECT_CI;
typedef MOL_PTR_VECT::iterator MOL_PTR_VECT_I;

}; // namespace RDKit

#endif
