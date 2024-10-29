#ifndef LAHUTA_RESIDUES_HPP
#define LAHUTA_RESIDUES_HPP

#include "GraphMol/Atom.h"
#include "GraphMol/RWMol.h"

namespace lahuta {

struct Residue {
  std::string chain_id;
  int number;
  std::string name;
  std::vector<const RDKit::Atom *> atoms;

  Residue() = default;
  Residue(const std::string &chain, int num, const std::string &name)
      : chain_id(chain), number(num), name(name) {}
};

class Residues {
public:
  explicit Residues(const RDKit::RWMol &mol_) : mol(mol_) { build_residues(mol_); }

  const std::vector<Residue> &get_residues() const { return residues_; }

  const std::unordered_map<std::string, std::vector<const Residue *>> &residue_map() const {
    return residues_by_name_;
  }

  /// filter residues by a predicate
  std::vector<Residue> filter(std::function<bool(const Residue &)> predicate) const;

  /// apply a function to each residue and collect the results
  template <typename ResultType>
  std::vector<ResultType> map(std::function<ResultType(const Residue &)> func) const;

  typedef std::vector<Residue>::const_iterator const_iterator;
  const_iterator begin() const { return residues_.begin(); }
  const_iterator end() const { return residues_.end(); }

private:
  void build_residues(const RDKit::RWMol &mol);

private:
  const RDKit::RWMol &mol;
  std::vector<Residue> residues_;
  std::unordered_map<std::string, std::vector<const Residue *>> residues_by_name_;
};

namespace residue_props {

template <typename ResultType> using RingProcFunc = std::function<ResultType(const Residue &, int)>;

inline Residue identity(const Residue &residue, int ring_size) { return residue; }

template <typename ResultType>
std::vector<ResultType>
get_aromatic_rings(const Residues &residues, RingProcFunc<ResultType> func = identity);

template <typename ReturnType>
ReturnType get_unknown_residues(const Residues &residues, const std::set<std::string> &KnownResiduesSet);

template <>
std::vector<Residue> get_unknown_residues<std::vector<Residue>>(
    const Residues &residues, const std::set<std::string> &KnownResiduesSet);

template <>
std::vector<int> get_unknown_residues<std::vector<int>>(
    const Residues &residues, const std::set<std::string> &KnownResiduesSet);

std::vector<std::vector<int>> tbl_find_aromatic_rings(const RDKit::RWMol &mol, const Residues &residues);

} // namespace residue_props

} // namespace lahuta

#endif // LAHUTA_RESIDUES_HPP
