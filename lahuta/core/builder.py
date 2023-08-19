"""Defines the LabeledNeighborPairsBuilder class.

Classes:
    LabeledNeighborPairsBuilder: A class to build LabeledNeighborPairs objects.
"""

from typing import Optional

import numpy as np
from Bio.Seq import Seq
from numpy.typing import NDArray

from lahuta.core.labeled_neighbors import LabeledNeighborPairs
from lahuta.lahuta_types.mdanalysis import AtomGroupType
from lahuta.msa.msa import MSAParser


class LabeledNeighborPairsBuilder:
    """Helper class to build LabeledNeighborPairs objects.

    The class provides a static method to map the pairs of atom indices to their corresponding atom names,
    residue IDs, and residue names. It also provides a static method to build a LabeledNeighborPairs object
    from the pairs of atom indices.

    Attributes:
        DTYPE (np.dtype): The data type of the LabeledNeighborPairs object.

    Static methods:
        map_pairs: Map the pairs of atom indices to their corresponding atom names, residue IDs, and residue names.
        build: Build a LabeledNeighborPairs object from the pairs of atom indices.

    """

    DTYPE = np.dtype(
        {
            "names": ["atom_names", "resids", "resnames"],
            "formats": ["<U10", "<U10", "<U10"],
        }
    )

    @staticmethod
    def map_pairs(atoms: AtomGroupType, seq: Optional[Seq] = None) -> tuple[NDArray[np.str_], ...]:
        """Map the pairs of atom indices to their corresponding atom names, residue IDs, and residue names.

        Args:
            atoms (AtomGroupType): The AtomGroup containing the atoms.
            seq (Optional[Seq], optional): The sequence of the AtomGroup. Defaults to None.

        Returns:
            tuple[NDArray[np.str_], ...]: A tuple containing the atom names, residue IDs, and residue names.
        """
        resindices = atoms.resindices
        if seq is None:
            resids = resindices.astype(np.str_)
        else:
            resids = np.array(MSAParser.to_indices_array(seq), dtype=np.str_)
            resids = resids[resindices]

        return atoms.names, atoms.resnames, resids

    @staticmethod
    def build(pairs: NDArray[np.int32], atoms: AtomGroupType, seq: Seq) -> "LabeledNeighborPairs":
        """Build a LabeledNeighborPairs object from the pairs of atom indices.

        Args:
            pairs (NDArray[np.int32]): The pairs of atom indices.
            atoms (AtomGroupType): The AtomGroup containing the atoms.
            seq (Seq): The sequence of the AtomGroup.

        Returns:
            LabeledNeighborPairs: A LabeledNeighborPairs object.
        """
        remapped = LabeledNeighborPairsBuilder.map_pairs(atoms, seq)

        names, resnames, resids = remapped
        data = np.empty(names.shape[0], dtype=LabeledNeighborPairsBuilder.DTYPE)
        data["atom_names"] = names
        data["resnames"] = resnames
        data["resids"] = resids

        return LabeledNeighborPairs(data[pairs])
