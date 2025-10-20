#ifndef LAHUTA_MODEL_RINGS_HPP
#define LAHUTA_MODEL_RINGS_HPP

#include <vector>

#include "constants.hpp"

namespace lahuta {

/// build aromatic ring vector for RingInfo
template <typename IndicesContainer> // IndicesContainer is an array with different sizes, hence the template
inline std::vector<int> make_aromatic_ring(int residue_start_idx, const IndicesContainer &indices) {
    std::vector<int> ring(indices.size());
    for (size_t k = 0; k < indices.size(); ++k) {
        ring[k] = residue_start_idx + indices[k];
    }
    return ring;
}

//
// Defs below were meant to be more performant replacement for how we handle rings,
// but they lead to weird lagging issues when run in parallel. At some point,
// we may want to revisit this.     - Besian, March 2025
//

constexpr int aroms_max_elem_idx      = 3; // the index of the last aromatic element in the array
constexpr int trp_aroms6_max_elem_idx = 4; // the index of the last aromatic element in the array

template <typename MoleculeType, typename RingContainer, typename BondContainer>
inline void add_phe_ring(MoleculeType& mol, int residue_start_idx,
                         RingContainer& aromatic_atom_indices,
                         BondContainer& aromatic_bond_indices,
                         int& max_arom_atom_idx, int& max_arom_bond_idx) {

    constexpr size_t ring_size = phe_arom_indices.size();
    std::array<int, ring_size> ring;
    std::array<int, ring_size> bonds;

    #pragma unroll
    for (size_t k = 0; k < ring_size; ++k) {
        ring[k] = residue_start_idx + phe_arom_indices[k];
    }

    #pragma unroll
    for (size_t ri = 0; ri < ring_size - 1; ++ri) {
        bonds[ri] = mol.addBond(ring[ri], ring[ri + 1], RDKit::Bond::AROMATIC);
    }
    bonds[ring_size - 1] = mol.addBond(ring[ring_size - 1], ring[0], RDKit::Bond::AROMATIC);

    aromatic_atom_indices.emplace_back(ring.begin(), ring.end());
    aromatic_bond_indices.emplace_back(bonds.begin(), bonds.end());

    max_arom_atom_idx = std::max(max_arom_atom_idx, ring[aroms_max_elem_idx]);
    max_arom_bond_idx = bonds[ring_size - 1];
}

template <typename MoleculeType, typename RingContainer, typename BondContainer>
inline void add_tyr_ring(MoleculeType& mol, int residue_start_idx, 
                         RingContainer& aromatic_atom_indices, 
                         BondContainer& aromatic_bond_indices,
                         int& max_arom_atom_idx, int& max_arom_bond_idx) {
    constexpr size_t ring_size = tyr_arom_indices.size();
    std::array<int, ring_size> ring;
    std::array<int, ring_size> bonds;

    #pragma unroll
    for (size_t k = 0; k < ring_size; ++k) {
        ring[k] = residue_start_idx + tyr_arom_indices[k];
    }

    #pragma unroll
    for (size_t ri = 0; ri < ring_size - 1; ++ri) {
        bonds[ri] = mol.addBond(ring[ri], ring[ri + 1], RDKit::Bond::AROMATIC);
    }
    bonds[ring_size - 1] = mol.addBond(ring[ring_size - 1], ring[0], RDKit::Bond::AROMATIC);

    aromatic_atom_indices.emplace_back(ring.begin(), ring.end());
    aromatic_bond_indices.emplace_back(bonds.begin(), bonds.end());

    max_arom_atom_idx = std::max(max_arom_atom_idx, ring[aroms_max_elem_idx]);
    max_arom_bond_idx = bonds[ring_size - 1];
}

template <typename MoleculeType, typename RingContainer, typename BondContainer>
inline void add_his_ring(MoleculeType& mol, int residue_start_idx, 
                         RingContainer& aromatic_atom_indices, 
                         BondContainer& aromatic_bond_indices,
                         int& max_arom_atom_idx, int& max_arom_bond_idx) {
    constexpr size_t ring_size = his_arom_indices.size();
    std::array<int, ring_size> ring;
    std::array<int, ring_size> bonds;

    #pragma unroll
    for (size_t k = 0; k < ring_size; ++k) {
        ring[k] = residue_start_idx + his_arom_indices[k];
    }

    #pragma unroll
    for (size_t ri = 0; ri < ring_size - 1; ++ri) {
        bonds[ri] = mol.addBond(ring[ri], ring[ri + 1], RDKit::Bond::AROMATIC);
    }
    bonds[ring_size - 1] = mol.addBond(ring[ring_size - 1], ring[0], RDKit::Bond::AROMATIC);


    aromatic_atom_indices.emplace_back(ring.begin(), ring.end());
    aromatic_bond_indices.emplace_back(bonds.begin(), bonds.end());

    max_arom_atom_idx = std::max(max_arom_atom_idx, ring[aroms_max_elem_idx]);
    max_arom_bond_idx = bonds[ring_size - 1];
}

template <typename MoleculeType, typename RingContainer, typename BondContainer>
inline void add_trp_rings(MoleculeType& mol, int residue_start_idx,
                         RingContainer& aromatic_atom_indices, 
                         BondContainer& aromatic_bond_indices,
                         int& max_arom_atom_idx, int& max_arom_bond_idx) {

    constexpr size_t ring_size5 = trp_arom_indices5.size();
    constexpr size_t ring_size6 = trp_arom_indices6.size();

    std::array<int, ring_size5> ring5;
    std::array<int, ring_size5> bonds5;

    std::array<int, ring_size6> ring6;
    std::array<int, ring_size6 - 1> bonds6; // skip the shared bond

    #pragma unroll
    for (size_t k = 0; k < ring_size5; ++k) {
        ring5[k] = residue_start_idx + trp_arom_indices5[k];
    }

    #pragma unroll
    for (size_t ri = 0; ri < ring_size5 - 1; ++ri) {
        bonds5[ri] = mol.addBond(ring5[ri], ring5[ri + 1], RDKit::Bond::AROMATIC);
    }
    bonds5[ring_size5 - 1] = mol.addBond(ring5[ring_size5 - 1], ring5[0], RDKit::Bond::AROMATIC);

    #pragma unroll
    for (size_t k = 0; k < ring_size6; ++k) {
        ring6[k] = residue_start_idx + trp_arom_indices6[k];
    }

    #pragma unroll
    for (size_t ri = 1; ri < ring_size6 - 1; ++ri) { // skip first bond (shared with the first ring)
        bonds6[ri-1] = mol.addBond(ring6[ri], ring6[ri + 1], RDKit::Bond::AROMATIC);
    }
    bonds6[ring_size6 - 2] = mol.addBond(ring6[ring_size6 - 1], ring6[0], RDKit::Bond::AROMATIC);

    aromatic_atom_indices.emplace_back(ring5.begin(), ring5.end());
    aromatic_atom_indices.emplace_back(ring6.begin() + 1, ring6.end());
    aromatic_bond_indices.emplace_back(bonds5.begin(), bonds5.end());
    aromatic_bond_indices.emplace_back(bonds6.begin(), bonds6.end());

    max_arom_atom_idx = std::max({max_arom_atom_idx, ring6[trp_aroms6_max_elem_idx]});
    max_arom_bond_idx = bonds6[ring_size6 - 2];
}

} // namespace lahuta

#endif // LAHUTA_MODEL_RINGS_HPP
