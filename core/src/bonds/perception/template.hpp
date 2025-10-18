#ifndef LAHUTA_BONDS_PERCEPTION_TEMPLATE_HPP
#define LAHUTA_BONDS_PERCEPTION_TEMPLATE_HPP

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Bond.h>
#include <rdkit/GraphMol/RWMol.h>

#include "key.hpp"
#include "span.hpp"

namespace lahuta::bonds {

// Order-independent residue atom signature used for consistency checks
struct AtomSignature {
  std::string name;
  unsigned int atomic_num = 0;

  bool operator==(const AtomSignature &other) const {
    return atomic_num == other.atomic_num && name == other.name;
  }

  bool operator!=(const AtomSignature &other) const { return !(*this == other); }
};

struct AtomTemplateData {
  int formal_charge = 0;
  bool is_aromatic = false;
};

struct BondTemplateData {
  unsigned int begin = 0;
  unsigned int end = 0;
  RDKit::Bond::BondType bond_type = RDKit::Bond::BondType::UNSPECIFIED;
  bool is_aromatic = false;
};

struct ResidueTemplate {
  std::vector<AtomTemplateData> atom_data;
  std::vector<BondTemplateData> bond_data;

  ResidueTemplate() = default;
  ResidueTemplate(std::vector<AtomTemplateData> atoms, std::vector<BondTemplateData> bonds)
      : atom_data(std::move(atoms)), bond_data(std::move(bonds)) {}
};

// One residue copy with indices into the source molecule
struct ResidueInstance {
  bonds::ResidueKey key;
  std::string instance_key;     // Cached string key for fast lookups
  std::vector<int> mol_indices; // Local indices within the molecule

  ResidueInstance() = default;

  explicit ResidueInstance(const bonds::ResidueKey &residue_key)
      : key(residue_key), instance_key(bonds::residue_key_utils::build_instance_key(residue_key)) {}
};

namespace detail {

// Create atom signatures for given atom indices
std::vector<AtomSignature> make_signature(const RDKit::RWMol &mol, span<const int> indices);

// Build a canonical, order-independent signature string
std::string make_signature_key(const RDKit::RWMol &mol, span<const int> indices);

// 64-bit order-independent signature hash
std::uint64_t make_signature_hash(const RDKit::RWMol &mol, span<const int> indices);

// Build a residue template from a small subgraph of source
std::optional<ResidueTemplate> build_residue_template(RDKit::RWMol &source, span<const int> indices);

// Apply template fields to the given instance (add/upgrade only)
bool apply_template_to_instance(
    RDKit::RWMol &target, const ResidueTemplate &templ, const ResidueInstance &instance);

} // namespace detail

} // namespace lahuta::bonds

#endif // LAHUTA_BONDS_PERCEPTION_TEMPLATE_HPP
