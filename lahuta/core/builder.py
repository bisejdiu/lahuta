"""Defines the LabeledNeighborPairsBuilder class.

Classes:
    LabeledNeighborPairsBuilder: A class to build LabeledNeighborPairs objects.
"""

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from numpy.typing import NDArray

from lahuta._types.mdanalysis import AtomGroupType
from lahuta.core.labeled_neighbors import LabeledNeighborPairs
from lahuta.msa.msa import MSAParser

__all__ = ["AtomMapper", "LabeledNeighborPairsBuilder"]


class AtomMapper:
    """Helper class to map atom indices to residue indices.

    The class provides a method to map atom indices to residue indices. It also provides a method to map
    residue indices to residue names and residue IDs.

    Attributes:
        atoms (AtomGroupType): An AtomGroupType object.
        prot (AtomGroupType): An AtomGroupType object containing only protein atoms.
        nonprot (AtomGroupType): An AtomGroupType object containing only non-protein atoms.

    Methods:
        map: Map atom indices to residue indices.
        factorize: Map residue indices to residue names and residue IDs.
    """

    def __init__(self, atoms: AtomGroupType):
        self.atoms = atoms
        self.prot, self.nonprot = self._mda_protein_select_split(atoms)

    @staticmethod
    def _mda_protein_select_split(atoms: AtomGroupType) -> tuple[AtomGroupType, AtomGroupType]:
        return atoms.select_atoms("protein"), atoms.select_atoms("not protein")

    def map(self, seq: Seq) -> NDArray[np.int32]:
        """Map atom indices to residue indices.

        Args:
            seq (Seq): A sequence as it results from the MSA.

        Returns:
            NDArray[np.int32]: Array of mapped residue indices.
        """
        prot_resindices = self._map_prot_resindices(seq)
        nonprot_resindices = self._map_nonprot_resindices(seq)

        mapped_resindices = self.sort_mapped_resindices(
            self.prot.resindices, self.nonprot.resindices, prot_resindices, nonprot_resindices
        )

        return mapped_resindices  # noqa: R504

    def _map_prot_resindices(self, seq: Seq) -> NDArray[np.int32]:
        mapped_prot_resindices = MSAParser.to_indices_array(seq)
        prot_resindices = self._factorize(self.prot.resindices)
        return mapped_prot_resindices[prot_resindices]

    def _map_nonprot_resindices(self, seq: Seq) -> NDArray[np.int32]:
        n_nonprot_residues = self.nonprot.residues.resindices.shape[0]
        nonprot_resindices = self._factorize(self.nonprot.resindices)
        shift_nonprot_resindices = np.arange(len(seq), len(seq) + n_nonprot_residues)
        return shift_nonprot_resindices[nonprot_resindices]

    @staticmethod
    def _factorize(resindices: NDArray[np.int32]) -> NDArray[np.int32]:
        return pd.factorize(resindices)[0]  # type: ignore

    def sort_mapped_resindices(
        self,
        prot_resindices: NDArray[np.int32],
        nonprot_resindices: NDArray[np.int32],
        mapped_prot_resindices: NDArray[np.int32],
        mapped_nonprot_resindices: NDArray[np.int32],
    ) -> NDArray[np.int32]:
        """Sort mapped residue indices.

        Args:
            prot_resindices (NDArray[np.int32]): Array of protein residue indices.
            nonprot_resindices (NDArray[np.int32]): Array of non-protein residue indices.
            mapped_prot_resindices (NDArray[np.int32]): Array of mapped protein residue indices.
            mapped_nonprot_resindices (NDArray[np.int32]): Array of mapped non-protein residue indices.

        Returns:
            NDArray[np.int32]: Sorted array of mapped residue indices.
        """
        indices = np.searchsorted(prot_resindices, nonprot_resindices)
        mapped_resindices = np.insert(mapped_prot_resindices, indices, mapped_nonprot_resindices)
        return mapped_resindices  # noqa: R504


class LabeledNeighborPairsBuilder:
    """Helper class to build LabeledNeighborPairs objects.

    The class provides a static method to map the pairs of atom indices to their corresponding atom names,
    residue IDs, and residue names. It also provides a static method to build a LabeledNeighborPairs object
    from the pairs of atom indices.

    Attributes:
        DTYPE (np.dtype): The data type of the LabeledNeighborPairs object.

    Methods:
        build: Build a LabeledNeighborPairs object from the pairs of atom indices.

    """

    DTYPE = np.dtype(
        {"names": ["chain_ids", "resids", "resnames", "names"], "formats": ["<U25", "<U25", "<U25", "<U25"]}
    )

    def __init__(self, atom_mapper: AtomMapper):
        self.atom_mapper = atom_mapper

    def build(self, pairs: NDArray[np.int32], seq: Seq) -> "LabeledNeighborPairs":
        """Build a LabeledNeighborPairs object from the pairs of atom indices and a sequence.

        Args:
            pairs (NDArray[np.int32]): Array of pairs of atom indices.
            seq (Seq): A sequence as it results from the MSA.

        Returns:
            LabeledNeighborPairs: A LabeledNeighborPairs object.
        """
        mapped_resindices = self.atom_mapper.map(seq)
        atoms = self.atom_mapper.atoms

        data = LabeledNeighborPairsBuilder.create_empty_struct_array(self.atom_mapper.atoms.n_atoms)
        data["chain_ids"] = atoms.chainIDs
        data["resnames"] = atoms.resnames
        data["resids"] = mapped_resindices
        data["names"] = atoms.names

        return LabeledNeighborPairs(data[pairs])

    @staticmethod
    def create_empty_struct_array(size: int) -> NDArray[np.void]:
        """Create an empty structured array.

        Args:
            size (int): The size of the array.

        Returns:
            NDArray[np.void]: An empty structured array.
        """
        return np.empty(size, dtype=LabeledNeighborPairsBuilder.DTYPE)
