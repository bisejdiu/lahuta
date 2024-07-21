#include <numeric>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>

#include "GraphMol/Atom.h"

// Custom hash function for a tuple of three strings
struct key_hash {
  std::size_t operator()(const std::tuple<std::string, std::string, std::string> &key) const {
    auto &[comp_id, atom_id1, atom_id2] = key;
    std::hash<std::string> hasher;
    return hasher(comp_id) ^ (hasher(atom_id1) << 1) ^ (hasher(atom_id2) << 2);
  }
};

template <typename T> std::set<T> union_of_many_sets(const std::initializer_list<std::set<T>> &sets) {
  std::set<T> result;
  for (const auto &set : sets) {
    result.insert(set.begin(), set.end());
  }
  return result;
}

// Define sets of amino acids, RNA bases, DNA bases, and peptide bases
inline std::set<std::string> amino_acid_names_long = {
    "HIS", "ARG", "LYS", "ILE", "PHE", "LEU",  "TRP",  "ALA", "MET", "PRO", "CYS", "ASN", "VAL", "GLY",
    "SER", "GLN", "TYR", "ASP", "GLU", "THR",  "SEC",  "PYL", "UNK", "MSE", "SEP", "TPO", "PTR", "PCA",
    "HYP", "HSD", "HSE", "HSP", "LSN", "ASPP", "GLUP", "HID", "HIE", "HIP", "LYN", "ASH", "GLH"};
inline std::set<std::string> amino_acid_names_dextro = {"DAL", "DAR"};
inline std::set<std::string> rna_base_names = {"A", "C", "T", "G", "I", "U", "N"};
inline std::set<std::string> dna_base_names = {"DA", "DC", "DT", "DG", "DI", "DU", "DN"};
inline std::set<std::string> peptide_base_names = {"APN", "CPN", "TPN", "GPN"};

inline std::set<std::string> combined_amino_acid_names = union_of_many_sets({amino_acid_names_long, amino_acid_names_dextro});
inline std::set<std::string> combined_base_names = union_of_many_sets({rna_base_names, dna_base_names, peptide_base_names});
inline std::set<std::string> combined_all_names = union_of_many_sets({combined_amino_acid_names, combined_base_names});

// Define the unordered_map using custom hash function
inline std::unordered_map<std::tuple<std::string, std::string, std::string>, int, key_hash> intra_bond_order_table;

// Initialize the table with precomputed keys
// FIX: turn into a constexpr function
inline void initialize_bond_order_table() {
  intra_bond_order_table[{"HIS", "CD2", "CG"}] = 2;
  intra_bond_order_table[{"HIS", "CE1", "ND1"}] = 2;
  intra_bond_order_table[{"ARG", "CZ", "NH2"}] = 2;
  intra_bond_order_table[{"PHE", "CE1", "CZ"}] = 2;
  intra_bond_order_table[{"PHE", "CD2", "CE2"}] = 2;
  intra_bond_order_table[{"PHE", "CD1", "CG"}] = 2;
  intra_bond_order_table[{"TRP", "CD1", "CG"}] = 2;
  intra_bond_order_table[{"TRP", "CD2", "CE2"}] = 2;
  intra_bond_order_table[{"TRP", "CE3", "CZ3"}] = 2;
  intra_bond_order_table[{"TRP", "CH2", "CZ2"}] = 2;
  intra_bond_order_table[{"ASN", "CG", "OD1"}] = 2;
  intra_bond_order_table[{"GLN", "CD", "OE1"}] = 2;
  intra_bond_order_table[{"TYR", "CD1", "CG"}] = 2;
  intra_bond_order_table[{"TYR", "CD2", "CE2"}] = 2;
  intra_bond_order_table[{"TYR", "CE1", "CZ"}] = 2;
  intra_bond_order_table[{"ASP", "CG", "OD1"}] = 2;
  intra_bond_order_table[{"GLU", "CD", "OE1"}] = 2;

  intra_bond_order_table[{"G", "C8", "N7"}] = 2;
  intra_bond_order_table[{"G", "C4", "C5"}] = 2;
  intra_bond_order_table[{"G", "C2", "N3"}] = 2;
  intra_bond_order_table[{"G", "C6", "O6"}] = 2;
  intra_bond_order_table[{"C", "C4", "N3"}] = 2;
  intra_bond_order_table[{"C", "C5", "C6"}] = 2;
  intra_bond_order_table[{"C", "C2", "O2"}] = 2;
  intra_bond_order_table[{"A", "C2", "N3"}] = 2;
  intra_bond_order_table[{"A", "C6", "N1"}] = 2;
  intra_bond_order_table[{"A", "C4", "C5"}] = 2;
  intra_bond_order_table[{"A", "C8", "N7"}] = 2;
  intra_bond_order_table[{"U", "C5", "C6"}] = 2;
  intra_bond_order_table[{"U", "C2", "O2"}] = 2;
  intra_bond_order_table[{"U", "C4", "O4"}] = 2;

  intra_bond_order_table[{"DG", "C8", "N7"}] = 2;
  intra_bond_order_table[{"DG", "C4", "C5"}] = 2;
  intra_bond_order_table[{"DG", "C2", "N3"}] = 2;
  intra_bond_order_table[{"DG", "C6", "O6"}] = 2;
  intra_bond_order_table[{"DC", "C4", "N3"}] = 2;
  intra_bond_order_table[{"DC", "C5", "C6"}] = 2;
  intra_bond_order_table[{"DC", "C2", "O2"}] = 2;
  intra_bond_order_table[{"DA", "C2", "N3"}] = 2;
  intra_bond_order_table[{"DA", "C6", "N1"}] = 2;
  intra_bond_order_table[{"DA", "C4", "C5"}] = 2;
  intra_bond_order_table[{"DA", "C8", "N7"}] = 2;
  intra_bond_order_table[{"DT", "C5", "C6"}] = 2;
  intra_bond_order_table[{"DT", "C2", "O2"}] = 2;
  intra_bond_order_table[{"DT", "C4", "O4"}] = 2;
}

int get_intra_bond_order(const std::string &comp_id, std::string atom_id1, std::string atom_id2);
RDKit::Atom::HybridizationType get_intra_bond_order(const RDKit::Atom *atom1, const RDKit::Atom *atom2);

