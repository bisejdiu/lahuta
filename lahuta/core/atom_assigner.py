from collections import OrderedDict

import numpy as np

from lahuta.config.atoms import PROT_ATOM_TYPES
from lahuta.config.smarts import ATOM_TYPES
from lahuta.core.assigners import (
    LegacyProteinTypeAssigner,
    VectorizedProteinTypeAssigner,
)
from lahuta.core.matchers import ParallelSmartsMatcher, SmartsMatcher


class AtomTypeAssigner:
    """
    A class for assigning atom types to atoms in a molecule.

    This class computes the atom types for a given molecule, utilizing various
    methods such as SMARTS pattern matching and protein atom type assignment.
    It can be configured to use different methods for SMARTS pattern matching
    (sequential or parallel) and protein atom type assignment (vectorized or legacy).
    """

    def __init__(self, mol, atomgroup, ta, parallel=False, legacy=True):
        self.mol = mol
        self.atomgroup = atomgroup
        print("->", self.atomgroup)
        self.protein_atomgroup = self.atomgroup.select_atoms("protein")
        self.ta = ta

        # self.atypes = OrderedDict(
        #     {x: i for i, x in enumerate(list(PROT_ATOM_TYPES.keys()))}
        # )
        self.atypes = OrderedDict(
            {
                "hbond acceptor": 0,
                "pos ionisable": 1,
                "carbonyl oxygen": 2,
                "weak hbond donor": 3,
                "carbonyl carbon": 4,
                "weak hbond acceptor": 5,
                "hbond donor": 6,
                "neg ionisable": 7,
                "aromatic": 8,
                "xbond acceptor": 9,
                "hydrophobe": 10,
            }
        )
        self.parallel = parallel
        self.legacy = legacy

        self.smarts_matcher_classes = {
            True: ParallelSmartsMatcher,
            False: SmartsMatcher,
        }

        self.protein_type_assigner_classes = {
            True: LegacyProteinTypeAssigner,
            False: VectorizedProteinTypeAssigner,
        }

    def _compute_smarts_types(self):
        """
        Compute atom types based on SMARTS pattern matching.

        Depending on the configuration, either the SmartsMatcher or the
        ParallelSmartsMatcher class will be used for SMARTS pattern matching.
        Returns an array of atom types as defined in the ATOM_TYPES dictionary.
        """
        smarts_matcher_class = self.smarts_matcher_classes[self.parallel]
        smarts_matcher = smarts_matcher_class(ATOM_TYPES)
        return smarts_matcher.compute(self.mol, self.atypes)

    def _compute_water_types(self, atypes_array):
        """
        Assign hydrogen bond donor and acceptor types to water molecules.

        Modifies the input atypes_array to set hydrogen bond donor and acceptor
        types for all water molecules in the atomgroup. Returns the modified
        atypes_array.
        """
        for atom in self.atomgroup.select_atoms(
            "resname SOL HOH TIP3 TIP4 WAT W and not name H*"
        ):
            atypes_array[atom.index, self.atypes["hbond acceptor"]] = 1
            atypes_array[atom.index, self.atypes["hbond donor"]] = 1

        return atypes_array

    def _compute_protein_types(self, atypes_array):
        """
        Assign protein atom types based on the chosen method.

        Depending on the configuration, either the VectorizedProteinTypeAssigner
        or the LegacyProteinTypeAssigner class will be used for protein atom
        type assignment. Returns an array of protein atom types as defined in
        the PROT_ATOM_TYPES dictionary.
        """
        protein_type_assigner_class = self.protein_type_assigner_classes[self.legacy]

        protein_type_assigner = protein_type_assigner_class(
            self.protein_atomgroup, self.ta
        )
        return protein_type_assigner.compute(atypes_array, self.atypes)

    def assign_atom_types(self):
        """
        Assign atom types to atoms in the molecule using the configured methods.

        Computes atom types for the molecule using the chosen methods for SMARTS
        pattern matching and protein atom type assignment. Returns an array of
        atom types for the entire molecule.
        """
        atypes_array = np.zeros((self.mol.NumAtoms(), len(PROT_ATOM_TYPES)))

        # atypes_array = self._compute_smarts_types()
        if self.atomgroup.n_atoms != self.protein_atomgroup.n_atoms:
            print("-> ", "not equal")
            atypes_array = self._compute_smarts_types()

        atypes_array = self._compute_water_types(atypes_array)
        atypes_array = self._compute_protein_types(atypes_array)

        return atypes_array
