#ifndef LAHUTA_CONVERT_HPP
#define LAHUTA_CONVERT_HPP

#include <gemmi/mmread_gz.hpp>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>
#include <vector>

using namespace gemmi;

namespace lahuta {

struct IR {
  std::vector<int> atom_indices;
  std::vector<int> atomic_numbers;
  std::vector<std::string> atom_names;
  std::vector<int> resids;
  std::vector<std::string> resnames;
  std::vector<std::string> chainlabels;
  std::vector<std::vector<float>> positions;

  IR() = delete;
  IR(const IR &) = default;
  IR(IR &&) = default;
  IR &operator=(const IR &) = default;
  IR &operator=(IR &&) = default;
  ~IR() = default;

  IR(std::vector<int> atom_indices, std::vector<int> atomic_numbers,
     std::vector<std::string> atom_names, std::vector<int> resids,
     std::vector<std::string> resnames, std::vector<std::string> chainlabels,
     std::vector<std::vector<float>> positions)
      : atom_indices(std::move(atom_indices)),
        atomic_numbers(std::move(atomic_numbers)),
        atom_names(std::move(atom_names)), resids(std::move(resids)),
        resnames(std::move(resnames)), chainlabels(std::move(chainlabels)),
        positions(std::move(positions)) {}
};

void IR_to_RWMol(RDKit::RWMol &mol, const IR &ir);
void gemmiStructureToRDKit(RDKit::RWMol &mol, const Structure &st,
                           RDKit::Conformer &conf, bool ign_h = true);

RDKit::RWMol filter_atoms(RDKit::RWMol &mol, std::vector<int> &indices);
RDKit::RWMol filter_with_conf(RDKit::RWMol &mol, std::vector<int> &indices);
RDKit::RWMol filter_with_bonds(const RDKit::RWMol &mol,
                               const std::vector<int> &indices);

} // namespace lahuta
//
#endif // LAHUTA_CONVERT_HPP
