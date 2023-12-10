"""Defines the LabeledNeighborPairsBuilder class.

Classes:
    LabeledNeighborPairsBuilder: A class to build LabeledNeighborPairs objects.
"""
import warnings
from typing import Optional

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from numpy.typing import NDArray

from lahuta._types.mdanalysis import AtomGroupType
from lahuta.core.neighbors import LabeledNeighborPairs
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
        self.sel_resindices = (self.atoms.universe.atoms - self.atoms).resindices

    @staticmethod
    def _mda_protein_select_split(atoms: AtomGroupType) -> tuple[AtomGroupType, AtomGroupType]:
        # TODO(bisejdiu): This selection won't work when we do not use synonyms when getting
        # the sequence from the Luni object. We need to synchronize these two methods.
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

        mapped_resindices = self.merge_mapped_indices(
            self.prot.resindices, self.nonprot.resindices, prot_resindices, nonprot_resindices
        )

        if self.sel_resindices.size > 0:
            mapped_resindices = self._mergemap_nonsel_resindices(mapped_resindices)

        return mapped_resindices

    def _map_prot_resindices(self, seq: Seq) -> NDArray[np.int32]:
        mapped_prot_resindices = MSAParser.to_indices_array(seq)
        prot_resindices = self._factorize(self.prot.resindices)
        return mapped_prot_resindices[prot_resindices]

    def _map_nonprot_resindices(self, seq: Seq) -> NDArray[np.int32]:
        n_nonprot_residues = self.nonprot.residues.resindices.shape[0]
        nonprot_resindices = self._factorize(self.nonprot.resindices)
        shift_nonprot_resindices = np.arange(len(seq), len(seq) + n_nonprot_residues)
        return shift_nonprot_resindices[nonprot_resindices]

    def _mergemap_nonsel_resindices(self, mapped_resindices: NDArray[np.int32]) -> NDArray[np.int32]:
        nonsel_indices = np.searchsorted(self.atoms.resindices, self.sel_resindices)
        with warnings.catch_warnings():
            # np.nan insertion is not supported by numpy.
            warnings.simplefilter("ignore")
            nonsel_resindices = np.full(nonsel_indices.shape, np.nan, dtype=float)
            mapped_resindices = np.insert(mapped_resindices, nonsel_indices, nonsel_resindices)

        return mapped_resindices

    @staticmethod
    def _factorize(resindices: NDArray[np.int32]) -> NDArray[np.int32]:
        return pd.factorize(resindices)[0]  # type: ignore

    def merge_mapped_indices(
        self,
        ref: NDArray[np.int32],
        target: NDArray[np.int32],
        mapped_ref: NDArray[np.int32],
        mapped_target: NDArray[np.int32],
    ) -> NDArray[np.int32]:
        """Merge two sets of mapped indices into a single, ordered array.

        Integrates two arrays of mapped indices, derived from reference (ref) and target arrays,
        The original ref and target arrays guide the merging process, as direct merging of the mapped arrays
        is not possible due to their transformed nature. The method ensures the final array
        respects the original order of the ref and target indices.

        Args:
            ref (NDArray[np.int32]): Array of reference indices.
            target (NDArray[np.int32]): Array of target indices.
            mapped_ref (NDArray[np.int32]): Mapped indices corresponding to the ref array.
            mapped_target (NDArray[np.int32]): Mapped indices corresponding to the target array.

        Returns:
            NDArray[np.int32]: A merged array containing elements from both mapped_ref and mapped_target in an order
                            that respects the original position of ref and target indices.

        The function identifies the appropriate insertion points for elements of the mapped_target array into the
        mapped_ref array, guided by the positions of target indices in relation to ref indices. This preserves the
        integrity of the original order in the combined mapped output.
        """
        insertion_points = np.searchsorted(ref, target)
        merged_indices = np.insert(mapped_ref, insertion_points, mapped_target)
        return merged_indices


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

    def build(self, pairs: NDArray[np.int32], seq: Seq, custom_fields: Optional[dict]=None) -> "LabeledNeighborPairs":
        """Build a LabeledNeighborPairs object from the pairs of atom indices, a sequence,
           and optional custom fields.

        Args:
            pairs (NDArray[np.int32]): Array of pairs of atom indices.
            seq (Seq): A sequence as it results from the MSA.
            custom_fields (dict, optional): A dictionary of custom fields and values.

        Returns:
            LabeledNeighborPairs: A LabeledNeighborPairs object.
        """
        mapped_resindices = self.atom_mapper.map(seq)
        atoms = self.atom_mapper.atoms.universe.atoms

        # Create the base structured array
        data = LabeledNeighborPairsBuilder.create_empty_struct_array(atoms.n_atoms)
        data["chain_ids"] = atoms.chainIDs
        data["resnames"] = atoms.resnames
        data["resids"] = mapped_resindices
        data["names"] = atoms.names

        # Extend the array if custom fields are provided
        if custom_fields:
            extended_dtype = self._extend_dtype_with_custom_fields(custom_fields)
            extended_data = np.empty(atoms.n_atoms, dtype=extended_dtype)

            # Copy data from the original array
            for field in LabeledNeighborPairsBuilder.DTYPE.names:
                extended_data[field] = data[field]

            # Add custom fields
            mapped_prot_resindices = MSAParser.to_indices_array(seq)
            for field_name, field_data in custom_fields.items():
                mapped_values = np.full(atoms.indices.shape, field_data.get("fill", ""), dtype="<U25")
                mapped_values[mapped_prot_resindices] = field_data["values"]
                extended_data[field_name] = mapped_values

            data = extended_data

        return LabeledNeighborPairs(data[pairs])

    @staticmethod
    def _extend_dtype_with_custom_fields(custom_fields: dict) -> np.dtype:
        """Extend the dtype with custom fields.

        Args:
            custom_fields (dict): A dictionary of custom fields and their values.

        Returns:
            np.dtype: An extended dtype.
        """
        new_fields = [(name, "<U25") for name in custom_fields]
        extended_dtype = np.dtype(list(LabeledNeighborPairsBuilder.DTYPE.descr) + new_fields)
        return extended_dtype

    @staticmethod
    def create_empty_struct_array(size: int) -> NDArray[np.void]:
        """Create an empty structured array.

        Args:
            size (int): The size of the array.

        Returns:
            NDArray[np.void]: An empty structured array.
        """
        return np.empty(size, dtype=LabeledNeighborPairsBuilder.DTYPE)
