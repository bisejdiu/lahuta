/*
    This code adapts functionality from the Mol* (mol-star) molecular visualization package,
    originally in TypeScript, for use in C++ with RDKit in the Lahuta project.

    Mol* is licensed under the MIT License (https://github.com/molstar/molstar),
    while this adaptation is released under the GNU General Public License (GPL).
*/
#ifndef LAHUTA_VALENCE_HPP
#define LAHUTA_VALENCE_HPP

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/ROMol.h>
#include <gemmi/elem.hpp>

using Bond = RDKit::Bond;
using HybridizationType = RDKit::Atom::HybridizationType;

namespace lahuta {

/**
 * Computes the formal charge, implicit hydrogen count, and hybridization type for each atom in a molecule.
 *
 * - Takes into account each atom's atomic number, degree, valence, and neighboring bonds.
 * - The behavior of the model can be configured via the `assign_charge` and `assign_h` options,
 *   which control whether formal charge and implicit hydrogens are automatically assigned.
 *
 * Usage:
 *   ValenceModel valence_model(true, true);
 *   valence_model.apply(mol);
 */
class ValenceModel {
public:
  ValenceModel(bool assign_charge = true, bool assign_h = true)
      : assign_charge_(assign_charge), assign_h_(assign_h) {}

  void apply(const RDKit::RWMol &mol) {
    for (auto &atom : mol.atoms()) {
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
   * Determines if the atom participates in a conjugated π-system, either through a double bond
   * or as N/O adjacent to a double bond, with the following exclusions:
   *
   * - N/O with degree 4 are not conjugated.
   * - N/O adjacent to P=O or S=O are excluded (e.g., sulfonamide N remains sp3).
   */
  bool is_conjugated(const RDKit::ROMol &mol, const RDKit::Atom &atom) const;

  // Returns the number of bonds the given atom has to atoms of the specified element.
  int get_element_count(RDKit::ROMol &mol, RDKit::Atom &atom, int element) const;

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

private:
  bool assign_charge_;
  bool assign_h_;
};

} // namespace lahuta

#endif // LAHUTA_VALENCE_HPP
