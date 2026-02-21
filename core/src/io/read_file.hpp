/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto to_s = [](auto&& arg) {
 *     if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, const char*>) return std::string(arg);
 *     else return std::string(arg);
 *   };
 *   return to_s("besian") + to_s("sejdiu") + to_s("@gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_READ_FILE_HPP
#define LAHUTA_READ_FILE_HPP

#include <gemmi/gz.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/mmread_gz.hpp>
#include <gemmi/pdb.hpp>
#include <rdkit/GraphMol/RWMol.h>

#include "io/convert.hpp"
#include "io/read_mmcif.hpp"

namespace lahuta {
using gemmi::CoorFormat;
namespace cif = gemmi::cif;

inline std::shared_ptr<RDKit::RWMol> make_molecule(cif::Document&& doc) {
  // mmCIF files for deposition may have more than one block:
  // coordinates in the first block and restraints in the others.
  for (size_t i = 1; i < doc.blocks.size(); ++i)
    if (doc.blocks[i].has_tag("_atom_site.id")) {
      gemmi::fail("2+ blocks are ok if only the first one has coordinates;\n _atom_site in block #" + std::to_string(i+1) + ": " + doc.source);
    }
  return make_mol_from_block(doc.blocks.at(0), doc.source);
}

template<typename T>
std::shared_ptr<RDKit::RWMol> read_and_make_molecule(T&& input) {
  auto format = gemmi::coor_format_from_ext(input.basepath());
  switch (format) {
    case CoorFormat::Pdb: {
      auto st = gemmi::read_pdb(input);
      return create_RDKit(st);
    }
    case CoorFormat::Mmcif:  return make_molecule(cif::read(input));
    case CoorFormat::Mmjson: return make_molecule(cif::read_mmjson(input));
    case CoorFormat::ChemComp: {
      auto st = gemmi::make_structure_from_chemcomp_doc(cif::read(input));
      return create_RDKit(st);
    }
    case CoorFormat::Unknown:
    case CoorFormat::Detect:
    default:
      gemmi::fail("Unsupported format for direct molecule creation: " + (input.path().empty() ? "coordinate file" : input.path()) + ".");
  }
  return nullptr;
}

} // namespace lahuta

#endif // LAHUTA_READ_FILE_HPP
