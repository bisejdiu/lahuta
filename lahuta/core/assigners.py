"""Contains the abstract base class (ABC) for assigning atom types to proteins.
It also contains two child classes that implement the abstract method from the ABC.

Classes:
    ```
    ProteinTypeAssignerBase: Abstract base class for assigning atom types to proteins.
    VectorizedProteinTypeAssigner: Efficient, vectorized assignment of atom types.
    LegacyProteinTypeAssigner: Traditional, loop-based assignment of atom types.
    ```

    
"""
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse import dok_matrix

from lahuta.config.atoms import PROT_ATOM_TYPES
from lahuta.config.smarts import AVAILABLE_ATOM_TYPES as ATypes
from lahuta.lahuta_types.mdanalysis import AtomGroupType

if TYPE_CHECKING:
    from numpy.typing import NDArray


class ProteinTypeAssignerBase(ABC):
    """Abstract Base Class for assigning atom types to proteins.

    This abstract base class (ABC) outlines the necessary structure and interface for child classes that handle
    the assignment of atom types to proteins.

    Attributes:
        protein_ag (AtomGroupType): Group of atoms in a protein that will be assigned atom types.

    Child classes:
        ```
        VectorizedProteinTypeAssigner: Efficient, vectorized assignment of atom types.
        LegacyProteinTypeAssigner: Traditional, loop-based assignment of atom types.
        ```
    """

    def __init__(self, protein_ag: AtomGroupType) -> None:
        self.protein_ag = protein_ag

    @abstractmethod
    def compute(self, atypes_array: dok_matrix) -> dok_matrix:
        """Abstract method to compute atom types.

        Must be implemented by child classes.

        Args:
            atypes_array (dok_matrix): Sparse array of atom types.

        Raises:
            NotImplementedError: If not implemented by child class.
        """
        raise NotImplementedError


class VectorizedProteinTypeAssigner(ProteinTypeAssignerBase):
    """Assigns atom types to proteins in a vectorized manner.

    Child class of ProteinTypeAssignerBase that uses NumPy array manipulations
    for efficient assignment of atom types.

    """

    def compute(self, atypes_array: dok_matrix) -> dok_matrix:
        """Compute atom types in a vectorized manner.

        Uses NumPy array manipulations for efficient assignment of atom types.

        Args:
            atypes_array (dok_matrix): Sparse array of atom types.

        Returns:
            (dok_matrix): Sparse array of assigned atom types.
        """
        resname_str = self.protein_ag.resnames.astype(str)
        atom_name_str = self.protein_ag.names.astype(str)
        atype_names = [member.lower() for member in list(ATypes)]

        atom_id_labels: NDArray[np.str_] = np.core.defchararray.add(
            np.core.defchararray.strip(resname_str),
            np.core.defchararray.strip(atom_name_str),
        )

        prot_atom_types_array = [list(PROT_ATOM_TYPES[key]) for key in atype_names]
        assert isinstance(atom_id_labels, np.ndarray)
        mask: NDArray[np.bool_] = np.array(
            [np.isin(atom_id_labels, prot_atom_types) for prot_atom_types in prot_atom_types_array]
        )

        true_indices = np.argwhere(mask)

        original_indices = self.protein_ag.indices[true_indices[:, 1]]
        atypes_array[original_indices, true_indices[:, 0]] = 1

        return atypes_array


class LegacyProteinTypeAssigner(ProteinTypeAssignerBase):
    """Assign atom types to proteins using a loop-based method.

    Child class of ProteinTypeAssignerBase that uses a traditional, loop-based approach for assignment of atom types.

    """

    def compute(self, atypes_array: dok_matrix) -> dok_matrix:
        """Compute atom types using a loop-based method.

        Uses a traditional, loop-based approach for assignment of atom types.

        Args:
            atypes_array (dok_matrix): Sparse array of atom types.

        Returns:
            (dok_matrix): Sparse array of assigned atom types.
        """
        for residue in self.protein_ag.residues:
            for atom in residue.atoms:
                for atom_type in list(PROT_ATOM_TYPES.keys()):
                    atypes_array[atom.index, ATypes[atom_type.upper()]] = 0

                for atom_type, atom_ids in PROT_ATOM_TYPES.items():
                    atom_id = residue.resname.strip() + atom.name.strip()
                    if atom_id in atom_ids:
                        atypes_array[atom.index, ATypes[atom_type.upper()]] = 1

        return atypes_array
