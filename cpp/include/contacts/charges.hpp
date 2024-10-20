/*
    This code adapts functionality from the Mol* (mol-star) molecular visualization package,
    originally in TypeScript, for use in C++ with RDKit in the Lahuta project.

    Mol* is licensed under the MIT License (https://github.com/molstar/molstar),
    while this adaptation is released under the GNU General Public License (GPL).
*/
#ifndef LAHUTA_CHARGES_HPP
#define LAHUTA_CHARGES_HPP

#include "atom_types.hpp"
#include "features.hpp"
#include <GraphMol/Atom.h>
#include <GraphMol/RWMol.h>

namespace lahuta {

using ResMap = std::unordered_map<std::string, std::vector<const RDKit::Atom *>>;

// FIX: Combine with "bonds/table.hpp" and put in one file
// FIX: These definitions are also not complete
const std::set<std::string> PositivelyChargedResidues = {"ARG", "HIS", "LYS"};
const std::set<std::string> PolymerNames = {"HIS", "ARG", "LYS", "ILE", "PHE", "LEU", "TRP",
                                            "ALA", "MET", "PRO", "CYS", "ASN", "VAL", "GLY",
                                            "SER", "GLN", "TYR", "ASP", "GLU", "THR"};
const std::set<std::string> ProteinBackboneAtoms = {"CA",  "C",   "N",   "O",   "O1",  "O2", "OC1",
                                                    "OC2", "OT1", "OT2", "OX1", "OXT", "H",  "H1",
                                                    "H2",  "H3",  "HA",  "HN",  "HXT", "BB"};

const std::set<std::string> NegativelyChargedResidues = {"GLU", "ASP"};
const std::set<std::string> BaseNames = {"DA", "DC", "DT", "DG", "DI", "DU", "DN"};

Feature create_feature(AtomType type, FeatureGroup group, const std::vector<const RDKit::Atom *> &members);

auto identify_feature_groups(const RDKit::RWMol &mol);
auto identify_negative_feature_groups(const RDKit::RWMol &mol);

[[nodiscard]] std::vector<Feature> add_positive_charges(const RDKit::RWMol &mol, ResMap &res_map);
[[nodiscard]] std::vector<Feature> add_negative_charges(const RDKit::RWMol &mol, ResMap &res_map);

} // namespace lahuta

#endif // LAHUTA_CHARGES_HPP
