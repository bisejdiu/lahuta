#include <algorithm>
#include <cstdint>

#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/MonomerInfo.h>

#include "bonds/bond_order.hpp"
#include "logging/logging.hpp"
#include "subgraph.hpp"
#include "template.hpp"

// clang-format off
namespace lahuta::bonds::detail {

std::vector<AtomSignature> make_signature(const RDKit::RWMol &mol, span<const int> indices) {
  std::vector<AtomSignature> signature;
  signature.reserve(indices.size());
  for (size_t j = 0; j < indices.size(); ++j) {
    int idx = indices[j];
    const RDKit::Atom *atom = mol.getAtomWithIdx(idx);
    const auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    std::string name = info ? residue_key_utils::trim_copy(info->getName()) : atom->getSymbol();
    signature.push_back(AtomSignature{name, static_cast<unsigned int>(atom->getAtomicNum())});
  }
  return signature;
}

std::string make_signature_key(const RDKit::RWMol &mol, span<const int> indices) {
  auto sig = make_signature(mol, indices);
  std::sort(sig.begin(), sig.end(), [](const AtomSignature &a, const AtomSignature &b) {
    if (a.name == b.name) return a.atomic_num < b.atomic_num;
    return a.name < b.name;
  });
  std::string key;
  key.reserve(sig.size() * 8);
  for (size_t i = 0; i < sig.size(); ++i) {
    if (i) key.push_back(';');
    key.append(sig[i].name);
    key.push_back('#');
    key.append(std::to_string(sig[i].atomic_num));
  }
  return key;
}

std::uint64_t make_signature_hash(const RDKit::RWMol &mol, span<const int> indices) {
  auto sig = make_signature(mol, indices);
  std::sort(sig.begin(), sig.end(), [](const AtomSignature &a, const AtomSignature &b) {
    if (a.name == b.name) return a.atomic_num < b.atomic_num;
    return a.name < b.name;
  });
  std::uint64_t hash = 1469598103934665603ull;
  auto fnv_step = [&](const char *data, size_t len) {
    for (size_t i = 0; i < len; ++i) {
      hash ^= static_cast<std::uint8_t>(data[i]);
      hash *= 1099511628211ull;
    }
  };
  for (const auto &a : sig) {
    fnv_step(a.name.data(), a.name.size());
    char buf[8];
    std::uint32_t z = static_cast<std::uint32_t>(a.atomic_num);
    buf[0] = static_cast<char> (z        & 0xFF);
    buf[1] = static_cast<char>((z >>  8) & 0xFF);
    buf[2] = static_cast<char>((z >> 16) & 0xFF);
    buf[3] = static_cast<char>((z >> 24) & 0xFF);
    fnv_step(buf, 4);
    buf[0] = 0xFF;
    fnv_step(buf, 1);
  }
  return hash;
}

std::optional<ResidueTemplate> build_residue_template(RDKit::RWMol &source, span<const int> indices) {
  if (indices.size() == 0) return std::nullopt;
  auto residue_mol = subgraph::build_rdkit_submol(source, indices, /*include_bonds=*/true);
  if (residue_mol.getNumAtoms() == 0) return std::nullopt;

  perceive_bond_orders_obabel(residue_mol);
  for (auto atom : residue_mol.atoms()) {
    atom->calcExplicitValence(false);
    if (atom->getAtomicNum() == 7 && atom->getFormalCharge() == 0 && atom->getExplicitValence() == 4) {
      atom->setFormalCharge(1);
    }
  }

  bool include_dative_bonds = true;
  RDKit::MolOps::symmetrizeSSSR(residue_mol, include_dative_bonds);
  RDKit::MolOps::setAromaticity(residue_mol);

  ResidueTemplate templ;
  templ.atom_data.reserve(residue_mol.getNumAtoms());
  for (const auto atom : residue_mol.atoms()) {
    AtomTemplateData d;
    d.formal_charge = atom->getFormalCharge();
    d.is_aromatic   = atom->getIsAromatic();
    templ.atom_data.push_back(d);
  }
  templ.bond_data.reserve(residue_mol.getNumBonds());
  for (const auto bond : residue_mol.bonds()) {
    BondTemplateData bd;
    bd.begin       = bond->getBeginAtomIdx();
    bd.end         = bond->getEndAtomIdx();
    bd.bond_type   = bond->getBondType();
    bd.is_aromatic = bond->getIsAromatic();
    templ.bond_data.push_back(bd);
  }
  return templ;
}

bool apply_template_to_instance(RDKit::RWMol &target, const ResidueTemplate &templ, const ResidueInstance &inst) {

  if (templ.atom_data.size() != inst.mol_indices.size()) {
    Logger::get_logger()->warn(
        "Residue {} {} has mismatched atom count (template={}, instance={})",
        inst.key.resname,
        inst.key.resnum,
        templ.atom_data.size(),
        inst.mol_indices.size());
    return false;
  }

  const auto n = templ.atom_data.size();

  // Apply atom-level properties
  for (size_t i = 0; i < n; ++i) {
    auto *atom = target.getAtomWithIdx(inst.mol_indices[i]);
    const auto &d = templ.atom_data[i];

    atom->setFormalCharge(d.formal_charge);
    atom->setIsAromatic  (d.is_aromatic);
  }

  // Add/update bonds
  for (const auto &bd : templ.bond_data) {
    if (bd.begin >= n || bd.end >= n) continue;

    int begin_idx = inst.mol_indices[bd.begin];
    int end_idx   = inst.mol_indices[bd.end];

    RDKit::Bond *bond = target.getBondBetweenAtoms(begin_idx, end_idx);
    if (!bond) {
      target.addBond(begin_idx, end_idx, bd.bond_type);
      bond = target.getBondBetweenAtoms(begin_idx, end_idx);
    } else {
      bond->setBondType(bd.bond_type);
    }

    bond->setIsAromatic(bd.is_aromatic);
  }

  return true;
}

} // namespace lahuta::bonds::detail
