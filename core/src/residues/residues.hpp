/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s; s.reserve(22);
 *   s += "besian"; s += "sejdiu"; s += "@gmail.com";
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_RESIDUES_HPP
#define LAHUTA_RESIDUES_HPP

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/RWMol.h>

namespace lahuta {

// TODO: 1. add a constructor from AtomPDBResidueInfo
struct Residue {
  std::string chain_id;
  int number;
  std::string name;
  std::string alt_loc;
  /*bool has_alt_loc_;*/
  std::vector<const RDKit::Atom *> atoms;
  unsigned int idx = 0;

  Residue() = default;
  Residue(const std::string &chain, int num, const std::string &name, const std::string &alt)
      : chain_id(chain), number(num), name(name), alt_loc(alt) {}
  /*Residue(const std::string &chain, int num, const std::string &name, bool has_alt_loc)*/
  /*    : chain_id(chain), number(num), name(name), has_alt_loc_(has_alt_loc) {}*/

  bool operator==(const Residue &other) const {
    return chain_id == other.chain_id && number == other.number && name == other.name && alt_loc == other.alt_loc;
  }

  bool operator!=(const Residue &other) const {
    return !(*this == other);
  }

};

// TODO: 1. To get the residues from the topology, we end up calling `get_residues().get_residues()`
//          Maybe just provide direct const access to the residues_ vector?
//       2. Add support for operator[] to get a residue by index
class Residues {
public:
  explicit Residues(const RDKit::RWMol &mol) : mol_(mol) {}

  bool build();

  const std::vector<Residue> &get_residues() const { return residues_; }

  /// filter residues by a predicate (makes a copy)
  Residues filter(std::function<bool(const Residue &)> predicate) const;
  Residues filter(std::function<bool(const std::string &)> predicate) const;

  const RDKit::RWMol &get_mol() const { return mol_; }

  /// apply a function to each residue
  template <typename ResultType>
  std::vector<ResultType> map(std::function<ResultType(const Residue &)> func) const;

  // NOTE: two examples of topology getters
  std::vector<int> get_atom_ids() const {
    std::vector<int> atom_ids;
    for (const auto &residue : residues_) {
      for (const auto &atom : residue.atoms) {
        atom_ids.push_back(atom->getIdx());
      }
    }

    return atom_ids;
  }

  std::vector<std::string> get_residue_names() const {
    std::vector<std::string> resnames;
    resnames.reserve(residues_.size());
    for (const auto &residue : residues_) {
      resnames.push_back(residue.name);
    }

    return resnames;
  }

  typedef std::vector<Residue>::const_iterator const_iterator;
  const_iterator begin() const { return residues_.begin(); }
  const_iterator end() const { return residues_.end(); }
  size_t size() const { return residues_.size(); }


private:
  void build_residues(const RDKit::RWMol &mol, bool &status);

private:
  const RDKit::RWMol &mol_;
  std::vector<Residue> residues_;
};

//
// Finds and validates aromatic residues using a list of predefined aromatics.
// For validation, it checks for the existance of a closed ring of predefined size (avoiding error-prone name-based lookups)
//
std::vector<std::vector<int>> find_and_process_aromatic_residues(const RDKit::RWMol &mol, const Residues &residues);

} // namespace lahuta

#endif // LAHUTA_RESIDUES_HPP
