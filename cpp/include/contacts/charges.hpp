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
#include "residues.hpp"
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

/*export const NucleicBackboneAtoms = new Set([*/
/*    'P', 'OP1', 'OP2', 'HOP2', 'HOP3',*/
/*    'O2\'', 'O3\'', 'O4\'', 'O5\'', 'C1\'', 'C2\'', 'C3\'', 'C4\'', 'C5\'',*/
/*    'H1\'', 'H2\'', 'H2\'\'', 'HO2\'', 'H3\'', 'H4\'', 'H5\'', 'H5\'\'', 'HO3\'', 'HO5\'',*/
/*    'O2*', 'O3*', 'O4*', 'O5*', 'C1*', 'C2*', 'C3*', 'C4*', 'C5*'*/
/*]);*/
const std::set<std::string> NucleicBackboneAtoms = {
    "P",  "OP1", "OP2", "HOP2", "HOP3", "O2\'", "O3\'", "O4\'", "O5\'", "C1\'", "C2\'",
    "C3\'", "C4\'", "C5\'", "H1\'", "H2\'", "H2\'\'", "HO2\'", "H3\'", "H4\'", "H5\'",
    "H5\'\'", "HO3\'", "HO5\'", "O2*", "O3*", "O4*", "O5*", "C1*", "C2*", "C3*", "C4*",
    "C5*"};

/*export const AminoAcidNamesL = new Set([*/
/*    'HIS', 'ARG', 'LYS', 'ILE', 'PHE', 'LEU', 'TRP', 'ALA', 'MET', 'PRO', 'CYS',*/
/*    'ASN', 'VAL', 'GLY', 'SER', 'GLN', 'TYR', 'ASP', 'GLU', 'THR', 'SEC', 'PYL',*/
/*    'UNK', // unknown amino acid from CCD*/
/*    'MSE', 'SEP', 'TPO', 'PTR', 'PCA', 'HYP', // common from CCD*/
/**/
/*    // charmm ff*/
/*    'HSD', 'HSE', 'HSP', 'LSN', 'ASPP', 'GLUP',*/
/**/
/*    // amber ff*/
/*    'HID', 'HIE', 'HIP', 'LYN', 'ASH', 'GLH',*/
/*]);*/

const std::set<std::string> AminoAcidNamesL = {
    "HIS", "ARG", "LYS", "ILE", "PHE", "LEU", "TRP", "ALA", "MET", "PRO", "CYS",
    "ASN", "VAL", "GLY", "SER", "GLN", "TYR", "ASP", "GLU", "THR", "SEC", "PYL",
    "UNK", // unknown amino acid from CCD
    "MSE", "SEP", "TPO", "PTR", "PCA", "HYP", // common from CCD

    // charmm ff
    "HSD", "HSE", "HSP", "LSN", "ASPP", "GLUP",

    // amber ff
    "HID", "HIE", "HIP", "LYN", "ASH", "GLH",
};

/*export const AminoAcidNamesD = new Set([*/
/*    'DAL', // D-ALANINE*/
/*    'DAR', // D-ARGININE*/
/*    'DSG', // D-ASPARAGINE*/
/*    'DAS', // D-ASPARTIC ACID*/
/*    'DCY', // D-CYSTEINE*/
/*    'DGL', // D-GLUTAMIC ACID*/
/*    'DGN', // D-GLUTAMINE*/
/*    'DHI', // D-HISTIDINE*/
/*    'DIL', // D-ISOLEUCINE*/
/*    'DLE', // D-LEUCINE*/
/*    'DLY', // D-LYSINE*/
/*    'MED', // D-METHIONINE*/
/*    'DPN', // D-PHENYLALANINE*/
/*    'DPR', // D-PROLINE*/
/*    'DSN', // D-SERINE*/
/*    'DTH', // D-THREONINE*/
/*    'DTR', // D-TRYPTOPHAN*/
/*    'DTY', // D-TYROSINE*/
/*    'DVA', // D-VALINE*/
/*    'DNE' // D-NORLEUCINE*/
/*    // ???  // D-SELENOCYSTEINE*/
/*]);*/
const std::set<std::string> AminoAcidNamesD = {
    "DAL", // D-ALANINE
    "DAR", // D-ARGININE
    "DSG", // D-ASPARAGINE
    "DAS", // D-ASPARTIC ACID
    "DCY", // D-CYSTEINE
    "DGL", // D-GLUTAMIC ACID
    "DGN", // D-GLUTAMINE
    "DHI", // D-HISTIDINE
    "DIL", // D-ISOLEUCINE
    "DLE", // D-LEUCINE
    "DLY", // D-LYSINE
    "MED", // D-METHIONINE
    "DPN", // D-PHENYLALANINE
    "DPR", // D-PROLINE
    "DSN", // D-SERINE
    "DTH", // D-THREONINE
    "DTR", // D-TRYPTOPHAN
    "DTY", // D-TYROSINE
    "DVA", // D-VALINE
    "DNE"  // D-NORLEUCINE
    // ???  // D-SELENOCYSTEINE
};

/*export const AminoAcidNames = SetUtils.unionMany(AminoAcidNamesL, AminoAcidNamesD);*/

const std::set<std::string> AminoAcidNames = {
    "HIS", "ARG", "LYS", "ILE", "PHE", "LEU", "TRP", "ALA", "MET", "PRO", "CYS",
    "ASN", "VAL", "GLY", "SER", "GLN", "TYR", "ASP", "GLU", "THR", "SEC", "PYL",
    "UNK", // unknown amino acid from CCD
    "MSE", "SEP", "TPO", "PTR", "PCA", "HYP", // common from CCD

    // charmm ff
    "HSD", "HSE", "HSP", "LSN", "ASPP", "GLUP",

    // amber ff
    "HID", "HIE", "HIP", "LYN", "ASH", "GLH",

    "DAL", // D-ALANINE
    "DAR", // D-ARGININE
    "DSG", // D-ASPARAGINE
    "DAS", // D-ASPARTIC ACID
    "DCY", // D-CYSTEINE
    "DGL", // D-GLUTAMIC ACID
    "DGN", // D-GLUTAMINE
    "DHI", // D-HISTIDINE
    "DIL", // D-ISOLEUCINE
    "DLE", // D-LEUCINE
    "DLY", // D-LYSINE
    "MED", // D-METHIONINE
    "DPN", // D-PHENYLALANINE
    "DPR", // D-PROLINE
    "DSN", // D-SERINE
    "DTH", // D-THREONINE
    "DTR", // D-TRYPTOPHAN
    "DTY", // D-TYROSINE
    "DVA", // D-VALINE
    "DNE"  // D-NORLEUCINE
    // ???  // D-SELENOCYSTEINE
};

auto identify_feature_groups(const RDKit::RWMol &mol);
auto identify_negative_feature_groups(const RDKit::RWMol &mol);

/*[[nodiscard]] std::vector<Feature> add_positive_charges(const RDKit::RWMol &mol, ResMap &res_map);*/
/*[[nodiscard]] std::vector<Feature> add_negative_charges(const RDKit::RWMol &mol, ResMap &res_map);*/
[[nodiscard]] FeatureVec add_positive_charges(const RDKit::RWMol &mol, const Residues &residues);
[[nodiscard]] FeatureVec add_negative_charges(const RDKit::RWMol &mol, const Residues &residues);

} // namespace lahuta

#endif // LAHUTA_CHARGES_HPP
