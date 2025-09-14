#ifndef LAHUTA_CONVERT_HPP
#define LAHUTA_CONVERT_HPP

#include <vector>

#include <gemmi/mmread_gz.hpp>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>

namespace lahuta {

struct IR {
  std::vector<int> atom_indices;
  std::vector<int> atomic_numbers;
  std::vector<std::string> atom_names;
  std::vector<int> resids;
  std::vector<std::string> resnames;
  std::vector<std::string> chainlabels;
  std::vector<std::vector<float>> positions;

  IR() = default;
  IR(const IR &) = default;
  IR(IR &&) = default;
  IR &operator=(const IR &) = default;
  IR &operator=(IR &&) = default;
  ~IR() = default;

  IR(std::vector<int> atom_indices, std::vector<int> atomic_numbers, std::vector<std::string> atom_names,
     std::vector<int> resids, std::vector<std::string> resnames, std::vector<std::string> chainlabels,
     std::vector<std::vector<float>> positions)
      : atom_indices(std::move(atom_indices)), atomic_numbers(std::move(atomic_numbers)),
        atom_names(std::move(atom_names)), resids(std::move(resids)), resnames(std::move(resnames)),
        chainlabels(std::move(chainlabels)), positions(std::move(positions)) {}
};

void IR_to_RWMol(RDKit::RWMol &mol, const IR &ir);
void add_atom_to_mol(
    RDKit::RWMol &mol, RDKit::Conformer &conf, const gemmi::Atom &atom, const gemmi::Residue &res,
    const gemmi::Chain &chain);
void create_RDKit_repr(RDKit::RWMol &mol, const gemmi::Structure &st, RDKit::Conformer &conf);
std::shared_ptr<RDKit::RWMol> create_RDKit(const gemmi::Structure &st);

} // namespace lahuta

#endif // LAHUTA_CONVERT_HPP
