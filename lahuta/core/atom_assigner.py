from typing import Dict, Type

import numpy as np
from numpy.typing import NDArray

from lahuta.config.atoms import PROT_ATOM_TYPES
from lahuta.config.smarts import AVAILABLE_ATOM_TYPES
from lahuta.core.assigners import (
    LegacyProteinTypeAssigner,
    VectorizedProteinTypeAssigner,
)
from lahuta.core.matchers import ParallelSmartsMatcher, SmartsMatcher, SmartsMatcherBase
from lahuta.types.mdanalysis import AtomGroupType
from lahuta.types.openbabel import MolType


class AtomTypeAssigner:
    """
    A class for assigning atom types to atoms in a molecule.

    This class computes the atom types for a given molecule, utilizing various
    methods such as SMARTS pattern matching and protein atom type assignment.
    It can be configured to use different methods for SMARTS pattern matching
    (sequential or parallel) and protein atom type assignment (vectorized or legacy).
    """

    def __init__(
        self,
        mda: AtomGroupType,
        mol: MolType,
        mapping: NDArray[np.int64],
        parallel: bool = False,
        legacy: bool = False,
    ) -> None:
        self.mda = mda
        self.mol = mol
        self.mapping = mapping
        self.protein_ag = self.mda.select_atoms("protein")

        self.atypes = AVAILABLE_ATOM_TYPES
        self.parallel = parallel
        self.legacy = legacy

        self.smarts_matcher_classes: Dict[bool, Type[SmartsMatcherBase]] = {
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
        Returns an array of atom types as defined in the SmartsPatternRegistry dictionary.
        """
        smarts_matcher_class = self.smarts_matcher_classes[self.parallel]
        smarts_matcher = smarts_matcher_class()
        return smarts_matcher.compute(self.mol)

    def _compute_water_types(self, atypes_array):
        """
        Assign hydrogen bond donor and acceptor types to water molecules.

        Modifies the input atypes_array to set hydrogen bond donor and acceptor
        types for all water molecules in the atomgroup. Returns the modified
        atypes_array.
        """
        water_ag = self.mda.select_atoms(
            "resname SOL HOH TIP3 TIP4 WAT W and not name H*"
        )
        # max_index = np.max(self.mda.universe.atoms.indices)
        # atom_mapping = np.full(max_index + 1, -1)
        # atom_mapping[self.mda.atoms.indices] = np.arange(self.mda.n_atoms)

        # TODO: vectorize this
        hbond_acceptor = self.atypes["hbond_acceptor".upper()]
        hbond_donor = self.atypes["hbond_donor".upper()]
        for atom in water_ag:
            atypes_array[self.mapping[atom.index], hbond_acceptor] = 1
            atypes_array[self.mapping[atom.index], hbond_donor] = 1

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

        protein_type_assigner = protein_type_assigner_class(self.protein_ag)
        return protein_type_assigner.compute(atypes_array)

    def assign_atom_types(self):
        """
        Assign atom types to atoms in the molecule using the configured methods.

        Computes atom types for the molecule using the chosen methods for SMARTS
        pattern matching and protein atom type assignment. Returns an array of
        atom types for the entire molecule.
        """
        atypes_array = np.zeros((self.mol.NumAtoms(), len(PROT_ATOM_TYPES)))

        # atypes_array = self._compute_smarts_types()
        if self.mda.n_atoms != self.protein_ag.n_atoms:
            atypes_array = self._compute_smarts_types()

        atypes_array = self._compute_water_types(atypes_array)
        atypes_array = self._compute_protein_types(atypes_array)

        return atypes_array
