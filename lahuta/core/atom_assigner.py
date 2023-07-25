from typing import Dict, Type

import numpy as np
from numpy.typing import NDArray

from lahuta.config.atoms import PROT_ATOM_TYPES
from lahuta.config.smarts import AVAILABLE_ATOM_TYPES
from lahuta.core.assigners import LegacyProteinTypeAssigner, VectorizedProteinTypeAssigner
from lahuta.core.matchers import ParallelSmartsMatcher, SmartsMatcher, SmartsMatcherBase
from lahuta.lahuta_types.mdanalysis import AtomGroupType
from lahuta.lahuta_types.openbabel import MolType


class AtomTypeAssigner:
    """
    Class for assigning atom types to atoms in a molecule.

    Handles the assignment of atom types to a given molecule.
    The class utilizes multiple methods such as SMARTS pattern matching and protein atom type assignment.
    It can be configured to use different methods for SMARTS pattern matching (sequential or parallel)
    and protein atom type assignment (vectorized or legacy).

    Attributes:
        mda (AtomGroupType): Atom group representing the molecular data.
        mol (MolType): The molecule to which atom types will be assigned.
        mapping (NDArray[np.int64]): Array for atom mapping.
        parallel (bool, optional): Flag to use parallel SMARTS pattern matching. Default is False.
        legacy (bool, optional): Flag to use legacy protein atom type assignment. Default is False.
        protein_ag (AtomGroupType): Atom group containing only protein atoms.
        atypes (dict): Available atom types.
        smarts_matcher_classes (dict): SMARTS pattern matching classes.
        protein_type_assigner_classes (dict): Protein atom type assigner classes.
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

    def _compute_smarts_types(self) -> NDArray[np.int8]:
        """
        Compute atom types based on SMARTS pattern matching.

        Depending on the configuration, either the SmartsMatcher or the ParallelSmartsMatcher
        class is used for SMARTS pattern matching.

        Returns:
            NDArray[np.int8]: Array of atom types as defined in the SmartsPatternRegistry dictionary.
        """
        smarts_matcher_class = self.smarts_matcher_classes[self.parallel]
        smarts_matcher = smarts_matcher_class()
        return smarts_matcher.compute(self.mol)

    def _compute_water_types(self, atypes_array: NDArray[np.int8]) -> NDArray[np.int8]:
        """
        Assign hydrogen bond donor and acceptor types to water molecules.

        This modifies the input atypes_array to set hydrogen bond donor and acceptor types
        for all water molecules in the atom group.

        Args:
            atypes_array (NDArray[np.int8]): Array of atom types.

        Returns:
            NDArray[np.int8]: Modified array of atom types with assigned types for water molecules.
        """
        water_ag = self.mda.select_atoms("resname SOL HOH TIP3 TIP4 WAT W and not name H*")
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

    def _compute_protein_types(self, atypes_array: NDArray[np.int8]) -> NDArray[np.int8]:
        """
        Assign protein atom types based on the chosen method.

        Depending on the configuration, either the VectorizedProteinTypeAssigner or the
        LegacyProteinTypeAssigner class is used for protein atom type assignment.

        Args:
            atypes_array (NDArray[np.int8]): Array of atom types.

        Returns:
            NDArray[np.int8]: Array of protein atom types as defined in PROT_ATOM_TYPES.
        """
        protein_type_assigner_class = self.protein_type_assigner_classes[self.legacy]

        protein_type_assigner = protein_type_assigner_class(self.protein_ag)
        return protein_type_assigner.compute(atypes_array)

    def assign_atom_types(self) -> NDArray[np.int8]:
        """
        Assign atom types to atoms in the molecule using the configured methods.

        Atom types for the molecule are computed using the chosen methods for SMARTS pattern matching and
        protein atom type assignment.

        Returns:
            NDArray[np.int8]: Array of atom types for the entire molecule.
        """
        atypes_array: NDArray[np.int8] = np.zeros((self.mol.NumAtoms(), len(PROT_ATOM_TYPES)), dtype=np.int8)

        # atypes_array = self._compute_smarts_types()
        if self.mda.n_atoms != self.protein_ag.n_atoms:
            atypes_array = self._compute_smarts_types()

        atypes_array = self._compute_water_types(atypes_array)
        atypes_array = self._compute_protein_types(atypes_array)

        return atypes_array
