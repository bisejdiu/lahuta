"""
Placeholder for the universe module.
"""

from typing import List, Tuple, Union, Dict
from typing_extensions import Protocol, Literal

import numpy as np
import pandas as pd
import MDAnalysis as mda
from openbabel import openbabel as ob

# from .protocol import ContactStrategy

from ..core.groups import AtomGroup
from ..core.universe import Universe
from ..core.neighbors import NeighborPairs
from ..utils.writers import DataFrameFactory

from ..utils.array import matching_array_indices


class CovalentContactStrategy:
    """A class to find covalent contacts between atoms in a molecule.

    Parameters
    ----------
    universe : mda.Universe or mda.AtomGroup
        The universe object or atom group to find contacts in.
    neighbor_pairs : NeighborPairs
        The neighbor pairs object.

    """

    def __init__(self, uniatom: Union[Universe, AtomGroup], neighbors: NeighborPairs):
        self._u = uniatom
        self.neighbors = neighbors

        self.pairs = self._pairs()
        self.distances = self._distances()

    @property
    def atoms(self) -> mda.AtomGroup:
        """Return an AtomGroup object for covalently bonded pairs.

        Returns
        -------
        ag : mda.AtomGroup
            An AtomGroup object containing the covalently bonded pairs.
        """
        return self._u.atoms[self.pairs]

    def is_covalent(self, pair: np.ndarray):
        """Check if two atoms are bonded.

        Parameters
        ----------
        pair : np.ndarray
            An array of shape (2,) where each row is a pair of atom indices.

        Returns
        -------
        boolean: bool
            True if the pair is bonded, False otherwise.

        """
        ix1, ix2 = pair
        for cov_bonded in ob.OBAtomAtomIter(
            self._u.atoms.universe.mol.GetAtom(int(ix1) + 1)
        ):
            if cov_bonded.GetId() == ix2:
                return True
        return False

    def contacts(
        self,
        output: Literal["dict", "dataframe"],
        dftype: Literal["compact", "expanded", "print"] = "compact",
    ) -> Union[Dict, pd.DataFrame]:
        """Get the covalent contacts.

        Parameters
        ----------
        output : Literal["dict", "dataframe"], optional
            The output of the output, by default "dict"

        Returns
        -------
        Union[Dict, pd.DataFrame]
            A dictionary/pandas DataFrame populated with covalent contacts.
        """

        # check input
        if output not in ["dict", "dataframe"]:
            raise ValueError(f"output must be 'dict' or 'dataframe', not {output}")

        if dftype not in ["compact", "expanded", "print"]:
            raise ValueError(
                f"dftype must be 'compact', 'expanded', or 'print', not {dftype}"
            )

        if output == "dict":
            df = DataFrameFactory(self.atoms, self.distances, dftype).dataframe()
            return df.to_dict(orient="list")
        else:
            return DataFrameFactory(self.atoms, self.distances, dftype).dataframe()

    def _pairs(self) -> np.ndarray:
        """Return the covalently bonded pairs."""
        cov_pair_indices = np.apply_along_axis(
            self.is_covalent, 1, self.neighbors.pairs
        )
        return self.neighbors.pairs[cov_pair_indices]

    def _distances(self):
        """Return distances between covalently bonded pairs."""
        matching_indices = matching_array_indices(self.neighbors.pairs, self.pairs)
        return self.neighbors.distances[matching_indices]

    @property
    def col1(self):
        """Get the first column of the covalent contact pairs.

        Returns
        -------
        col1 : np.ndarray
            An array of shape (n_pairs,) containing the atom indices of the first column of covalently bonded pairs.
        """
        return self.atoms[:, 0]

    @property
    def col2(self):
        """Get the second column of the covalent contact pairs.

        Returns
        -------
        col2 : np.ndarray
            An array of shape (n_pairs,) containing the atom indices of the second column of covalently bonded pairs.
        """
        return self.atoms[:, 1]
