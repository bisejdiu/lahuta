#ifndef LAHUTA_READ_MMCIF_HPP
#define LAHUTA_READ_MMCIF_HPP

#include <stdexcept>
#include <vector>

#include <gemmi/atox.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/elem.hpp>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>

#include "gemmi/numb.hpp"
#include "logging.hpp"

// clang-format off
namespace lahuta {

inline std::shared_ptr<RDKit::RWMol> make_mol_from_block(const gemmi::cif::Block &block_) {
  namespace cif = gemmi::cif;

  auto &block = const_cast<cif::Block &>(block_);
  auto mol    = std::make_shared<RDKit::RWMol>();

  enum {
    kId=0, kGroupPdb, kSymbol, kLabelAtomId, kAltId, kLabelCompId,
    kLabelAsymId, kLabelEntityId, kLabelSeqId, kInsCode,
    kX, kY, kZ, kOcc, kBiso, kCharge,
    kAuthSeqId, kAuthCompId, kAuthAsymId, kAuthAtomId, kModelNum
  };

  cif::Table atom_table = block.find("_atom_site.",
                                {"id",
                                 "?group_PDB",
                                 "type_symbol",
                                 "?label_atom_id",
                                 "label_alt_id",
                                 "?label_comp_id",
                                 "label_asym_id",
                                 "?label_entity_id",
                                 "?label_seq_id",
                                 "?pdbx_PDB_ins_code",
                                 "Cartn_x",
                                 "Cartn_y",
                                 "Cartn_z",
                                 "occupancy",
                                 "B_iso_or_equiv",
                                 "?pdbx_formal_charge",
                                 "auth_seq_id",
                                 "?auth_comp_id",
                                 "?auth_asym_id",
                                 "?auth_atom_id",
                                 "?pdbx_PDB_model_num"
                                });

  if (atom_table.length() == 0) throw std::runtime_error("No _atom_site table found in the mmCIF file");

  const int kAsymId = atom_table.first_of(kAuthAsymId, kLabelAsymId);
  const int kCompId = atom_table.first_of(kAuthCompId, kLabelCompId);
  const int kAtomId = atom_table.first_of(kAuthAtomId, kLabelAtomId);

  if (!atom_table.has_column(kCompId)) throw std::runtime_error("Neither _atom_site.label_comp_id nor auth_comp_id found");
  if (!atom_table.has_column(kAtomId)) throw std::runtime_error("Neither _atom_site.label_atom_id nor auth_atom_id found");

  // we expect that models (if multiple) are ordered by model number (1, 1, 2, 2, 3, 3)
  std::vector<std::vector<RDGeom::Point3D>> model_coords;
  model_coords.reserve(20);

  std::string current_model;
  for (auto row : atom_table) {
    std::string model_id = row.has(kModelNum) ? row.str(kModelNum) : "1";

    if (model_coords.empty() || model_id != current_model) {
      current_model = model_id;
      model_coords.emplace_back();
      model_coords.back().reserve(atom_table.length());
    }

    model_coords.back().emplace_back(cif::as_number(row[kX]), cif::as_number(row[kY]), cif::as_number(row[kZ]));

    if (model_coords.size() > 1) continue;

    std::string symbol = cif::as_string(row[kSymbol]);
    int atomic_number  = gemmi::Element(symbol).atomic_number();
    auto *rdkit_atom   = new RDKit::Atom(atomic_number);

    if (row.has2(kCharge)) rdkit_atom->setFormalCharge(cif::as_int(row[kCharge]));

    mol->addAtom(rdkit_atom, true, true);

    auto atom_name   = cif::as_string(row[kAtomId]);
    auto res_name    = cif::as_string(row[kCompId]);
    auto chain_name  = cif::as_string(row[kAsymId]);
    auto auth_seq_id = cif::as_string(row[kAuthSeqId]);
    int  seq_id      = gemmi::string_to_int(auth_seq_id, false);
    char alt_loc     = cif::as_char(row[kAltId], '\0');
    auto alt_loc_str = (alt_loc == '\0') ? "" : std::string(1, alt_loc);
    int  serial      = gemmi::string_to_int(row[kId], false);

    auto *info = new RDKit::AtomPDBResidueInfo(atom_name, serial, alt_loc_str, res_name, seq_id, chain_name);

    if (row.has2(kGroupPdb)) {
      char het_flag = '\0';
      for (int i = 0; i < 2; ++i) {
        const char c = toupper(row[kGroupPdb][i]);
        if (c == 'A' || c == 'H') het_flag = c;
      }
      info->setIsHeteroAtom(het_flag == 'H'); // it is false by default
    }

    if (row.has2(kLabelSeqId)) {
      int label_seq = cif::as_int(row[kLabelSeqId]);
      info->setResidueNumber(label_seq);
    }

    info->setMonomerType(RDKit::AtomMonomerInfo::PDBRESIDUE);
    rdkit_atom->setMonomerInfo(info);
  }

  const std::size_t atoms_per_model = model_coords[0].size();

  auto new_end = std::remove_if(model_coords.begin(), model_coords.end(), [atoms_per_model](auto const &coords) {
      return coords.size() != atoms_per_model;
    }
  );

  if (new_end != model_coords.end()) {
    model_coords.erase(new_end, model_coords.end());
    Logger::get_logger()->warn("Dropping {} models with inconsistent number of atoms", model_coords.size());
  }

  unsigned int conf_id = 0;
  for (auto& coords : model_coords) {
    RDKit::Conformer *conformer = new RDKit::Conformer(atoms_per_model);
    conformer->setId(conf_id++);
    conformer->setAllAtomPositions(std::move(coords));
    mol->addConformer(conformer, true);
  }

  if (mol->getNumConformers() != 0) {
    Logger::get_logger()->debug("Loaded {} models with {} atoms", mol->getNumConformers(), mol->getNumAtoms());
  }

  mol->updatePropertyCache(false);
  return mol;
}

} // namespace lahuta

#endif // LAHUTA_READ_MMCIF_HPP
