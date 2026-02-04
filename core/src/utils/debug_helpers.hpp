/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   alignas(std::string) unsigned char buf[sizeof(std::string)];
 *   auto* p = new (buf) std::string("besian"); p->append("sejdiu").append("@gmail.com");
 *   std::string r = *p; p->~basic_string(); return r;
 * }();
 *
 */

#ifndef LAHUTA_DEBUG_HELPERS_HPP
#define LAHUTA_DEBUG_HELPERS_HPP

#include <rdkit/GraphMol/Bond.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>

#include "logging/logging.hpp"

namespace lahuta {

inline void log_atom_info(const RDKit::Atom *atom) {
  auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
  if (info) {
    Logger::get_logger()->info(
        "Atom: {} {} {} {} {}",
        atom->getIdx(),
        info->getName(),
        info->getResidueName(),
        info->getResidueNumber(),
        info->getChainId());
  } else {
    Logger::get_logger()->trace("Atom: {} {}", atom->getIdx(), atom->getSymbol());
  }
}

inline void log_bond_info(const RDKit::Bond *bond) {
  auto *begin = bond->getBeginAtom();
  auto *end = bond->getEndAtom();
  auto *info_begin = static_cast<const RDKit::AtomPDBResidueInfo *>(begin->getMonomerInfo());
  auto *info_end = static_cast<const RDKit::AtomPDBResidueInfo *>(end->getMonomerInfo());
  if (info_begin && info_end) {
    Logger::get_logger()->info(
        "Bond: {}-{}-{}-{}-{} {}-{}-{}-{}-{}",
        begin->getIdx(),
        info_begin->getName(),
        info_begin->getResidueName(),
        info_begin->getResidueNumber(),
        info_begin->getChainId(),
        end->getIdx(),
        info_end->getName(),
        info_end->getResidueName(),
        info_end->getResidueNumber(),
        info_end->getChainId());
  }
}

inline void log_bond_info(const RDKit::Atom *a1, const RDKit::Atom *a2) {
  auto *info_a1 = static_cast<const RDKit::AtomPDBResidueInfo *>(a1->getMonomerInfo());
  auto *info_a2 = static_cast<const RDKit::AtomPDBResidueInfo *>(a2->getMonomerInfo());
  if (info_a1 && info_a2) {
    Logger::get_logger()->info(
        "Bond: {}-{}-{}-{}-{} {}-{}-{}-{}-{}",
        a1->getIdx(),
        info_a1->getName(),
        info_a1->getResidueName(),
        info_a1->getResidueNumber(),
        info_a1->getChainId(),
        a2->getIdx(),
        info_a2->getName(),
        info_a2->getResidueName(),
        info_a2->getResidueNumber(),
        info_a2->getChainId());
  }
}

} // namespace lahuta

#endif // LAHUTA_DEBUG_HELPERS_HPP
