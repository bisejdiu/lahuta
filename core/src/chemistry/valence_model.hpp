#ifndef LAHUTA_VALENCE_HPP
#define LAHUTA_VALENCE_HPP

#include <gemmi/elem.hpp>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Bond.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/PeriodicTable.h>
#include <rdkit/GraphMol/ROMol.h>

using Bond              = RDKit::Bond;
using HybridizationType = RDKit::Atom::HybridizationType;

namespace lahuta {

/**
 * Computes the formal charge, implicit hydrogen count, and hybridization type for each atom in a molecule.
 *
 * - Takes into account each atom's atomic number, degree, valence, and neighboring bonds.
 *
 * Usage:
 *   ValenceModel valence_model;
 *   valence_model.apply(mol);
 */
class ValenceModel {
public:
  void apply(const RDKit::RWMol &mol) {
    for (const auto atom : mol.atoms()) {
      molstar_valence_model(mol, *atom);
    }
  }

  /// Returns true if the given atom has a double or aromatic bond.
  bool has_double_or_aromatic_bond(const RDKit::Atom &atom, const RDKit::ROMol &mol) const;

  /// exclude N/O adjacent to P=O or S=O
  bool is_excluded_bond(const RDKit::Atom &atom_b, const RDKit::Atom &atom_c) const;

  /// True if the atom has a neighbor with a double or aromatic bond that is not excluded.
  bool is_conjugated_and_not_excluded(const RDKit::Atom &atom_b, const RDKit::ROMol &mol) const;

  /**
   * Determines if the atom participates in a conjugated pi-system, either through a double bond
   * or as N/O adjacent to a double bond, with the following exclusions:
   *
   * - N/O with degree 4 are not conjugated.
   * - N/O adjacent to P=O or S=O are excluded (e.g., sulfonamide N remains sp3).
   */
  bool is_conjugated(const RDKit::ROMol &mol, const RDKit::Atom &atom) const;

  // Returns the number of bonds the given atom has to atoms of the specified element.
  int get_element_count(const RDKit::ROMol &mol, const RDKit::Atom &atom, int element) const;

  /// True if the atom is bound to sulfur or a metal.
  bool is_bound_to_sulfur_or_metal(const RDKit::ROMol &mol, RDKit::Atom *atom) const;

  /// True if the atom is an amidine or guanidine nitrogen.
  bool is_amidine_or_guanidine_nitrogen(int degree, int h_count, int valence) const;

  /// True if the atom has a neighbor with a double-bonded oxygen.
  bool has_neighbor_with_double_bonded_oxygen(const RDKit::ROMol &mol, RDKit::Atom *atom) const;

  /// Assigns the hybridization type based on the total coordination.
  HybridizationType assign_geometry(int total_coordination) const;

  /// Assigns the formal charge, implicit hydrogen count, and hybridization type for the given atom.
  void molstar_valence_model(const RDKit::ROMol &mol, RDKit::Atom &atom);
};

} // namespace lahuta

#endif // LAHUTA_VALENCE_HPP
