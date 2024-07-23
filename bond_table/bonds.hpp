#include <set>
#include <string>

#include "bond_order_hash.hpp"
#include "res_name_table.hpp"
#include "tokens_table.cpp"

#include "GraphMol/Atom.h"

template <typename T>
std::set<T> union_of_many_sets(const std::initializer_list<std::set<T>> &sets) {
  std::set<T> result;
  for (const auto &set : sets) {
    result.insert(set.begin(), set.end());
  }
  return result;
}

inline std::set<std::string> combined_base_names = {
    "A",  "C",  "T",  "G",  "I",  "U",   "N",   "DA",  "DC",
    "DT", "DG", "DI", "DU", "DN", "APN", "CPN", "TPN", "GPN"};

int get_intra_bond_order(const std::string &comp_id,
                         const std::string *atom_id1,
                         const std::string *atom_id2);

RDKit::Atom::HybridizationType get_intra_bond_order(const RDKit::Atom *atom1,
                                                    const RDKit::Atom *atom2);

inline auto get_bond_order(const std::string &comp_id,
                          const std::string &atom_id1,
                          const std::string &atom_id2) {
  std::string key = comp_id + "_" + atom_id1 + "_" + atom_id2;
  auto entry = BondOrderTable::get_bond_order(key.c_str(), key.length());
  return entry;
}

inline auto getResName(const std::string &comp_id) {
  auto entry = ResNameTable::getResName(comp_id.c_str(), comp_id.length());
  return entry;
}

inline auto getToken(const std::string &comp_id) {
  auto entry = TokenTable::getToken(comp_id.c_str(), comp_id.length());
  return entry;
}
