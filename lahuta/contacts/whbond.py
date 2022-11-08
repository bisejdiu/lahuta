"""
Placeholder for the universe module.
"""

from typing import Union, Dict
from typing_extensions import Literal

import numpy as np
import pandas as pd

from ..core.groups import AtomGroup
from ..core.universe import Universe
from ..core.neighbors import NeighborPairs, HBondNeighborPairs, NeighborPairsBase
from ..utils.writers import DataFrameFactory
from ..utils.array import matching_array_indices

from ..config import config


class WeakHBondContactStrategy:
    """A class to find hbond contacts between atoms in a molecule.

    Parameters
    ----------
    universe : mda.Universe
        The universe object.
    neighbor_pairs : NeighborPairs
        The neighbor pairs object.

    """

    vdw_comp_factor = 0.1

    def __init__(self, uniatom: Universe, neighbors: NeighborPairsBase):
        self._u = uniatom
        self.neighbors = neighbors

        self.pairs = self._pairs()
        self.distances = self._distances()

    @property
    def atoms(self) -> AtomGroup:
        """Return an AtomGroup object for metal contacts.

        Returns
        -------
        ag : AtomGroup
            An AtomGroup object containing the metal contacts.
        """
        return self._u.atoms[self.pairs]

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

    def _pairs(self):
        """Get the metal contacts.

        Returns
        -------
        pairs : np.ndarray
            An array of shape (n_pairs, 2) where each row is a pair of atom indices.
        """

        hbond_atom12 = (
            self.neighbors.type_filter("hbond acceptor", 0)
            .type_filter("weak hbond donor", 1)
            .contact_type("hbond")
            .hbond_distance_filter(col=0)
            .hbond_angle_filter(col=1, weak=True)
        )

        hbond_atom21 = (
            self.neighbors.type_filter("hbond acceptor", 1)
            .type_filter("weak hbond donor", 0)
            .contact_type("hbond")
            .hbond_distance_filter(col=1)
            .hbond_angle_filter(col=0, weak=True)
        )

        return np.concatenate((hbond_atom12.pairs, hbond_atom21.pairs), axis=0)

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
