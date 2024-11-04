/*
    This code adapts functionality from the Mol* (mol-star) molecular visualization package,
    originally in TypeScript, for use in C++ with RDKit in the Lahuta project.

    Mol* is licensed under the MIT License (https://github.com/molstar/molstar),
    while this adaptation is released under the GNU General Public License (GPL).
*/
#ifndef LAHUTA_CHARGES_HPP
#define LAHUTA_CHARGES_HPP

#include "entities.hpp"
#include "residues.hpp"
#include <GraphMol/Atom.h>
#include <GraphMol/RWMol.h>

namespace lahuta {

// FIX: Combine with "bonds/table.hpp" and put in one file
// FIX: These definitions are also not complete

const std::set<std::string> AminoAcidNamesL = {
    "HIS",
    "ARG",
    "LYS",
    "ILE",
    "PHE",
    "LEU",
    "TRP",
    "ALA",
    "MET",
    "PRO",
    "CYS",
    "ASN",
    "VAL",
    "GLY",
    "SER",
    "GLN",
    "TYR",
    "ASP",
    "GLU",
    "THR",
    "SEC",
    "PYL",
    "UNK", // unknown amino acid from CCD
    "MSE",
    "SEP",
    "TPO",
    "PTR",
    "PCA",
    "HYP", // common from CCD

    // charmm ff
    "HSD",
    "HSE",
    "HSP",
    "LSN",
    "ASPP",
    "GLUP",

    // amber ff
    "HID",
    "HIE",
    "HIP",
    "LYN",
    "ASH",
    "GLH",
};

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
};

const std::set<std::string> AminoAcidNames = {
    "HIS",
    "ARG",
    "LYS",
    "ILE",
    "PHE",
    "LEU",
    "TRP",
    "ALA",
    "MET",
    "PRO",
    "CYS",
    "ASN",
    "VAL",
    "GLY",
    "SER",
    "GLN",
    "TYR",
    "ASP",
    "GLU",
    "THR",
    "SEC",
    "PYL",
    "UNK", // unknown amino acid from CCD
    "MSE",
    "SEP",
    "TPO",
    "PTR",
    "PCA",
    "HYP", // common from CCD

    // charmm ff
    "HSD",
    "HSE",
    "HSP",
    "LSN",
    "ASPP",
    "GLUP",

    // amber ff
    "HID",
    "HIE",
    "HIP",
    "LYN",
    "ASH",
    "GLH",

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
};

auto identify_positive_charge_groups(const RDKit::RWMol &mol);
auto identify_negative_charge_groups(const RDKit::RWMol &mol);

[[nodiscard]] GroupEntityCollection add_positive_charges(const RDKit::RWMol &mol, const Residues &residues);
[[nodiscard]] GroupEntityCollection add_negative_charges(const RDKit::RWMol &mol, const Residues &residues);

} // namespace lahuta

#endif // LAHUTA_CHARGES_HPP
