"""
Placeholder
"""

from abc import ABC, abstractmethod
from typing import Any, Dict, Literal, Union

import numpy as np
import pandas as pd

from ..core.groups import AtomGroup
from ..core.neighbors import NeighborPairs
from ..core.universe import Universe
from ..utils.array_utils import matching_indices
from ..utils.writers import DataFrameFactory

# class ContactStrategy(Protocol):
#     """A protocol for contact strategies."""

#     universe: mda.Universe
#     neighbor_pairs: NeighborPairs

#     @property
#     def contact_pairs(self, **kwargs) -> np.ndarray:
#         """Retrieve the contact pairs indices."""

#     @property
#     def contact_atoms(self, **kwargs) -> mda.AtomGroup:
#         """Retrieve the contact atoms."""

#     def _contact_pairs(self, **kwargs) -> np.ndarray:
#         """Main method for computing the contact pairs."""


class ContactBase(ABC):
    """Abstract method for computing contacts."""

    def __init__(
        self,
        ua: Union[Universe, AtomGroup],
        neighbor_pairs: NeighborPairs,
        **kwargs,
    ):
        self._ua = ua
        self.neighbors = neighbor_pairs

        self.pairs = self.compute_contacts(**kwargs)
        self.distances = self.filter_distances()

    @property
    def indices(self, **kwargs) -> np.ndarray:
        """Retrieve the contact pairs indices."""
        return np.sort(np.unique([self.pairs[:, 0], self.pairs[:, 1]]))

    @property
    def atoms(self) -> AtomGroup:
        """Retrieve the contact atoms."""
        return self._ua.atoms[self.pairs]

    @property
    def col1(self):
        """Return the first column of the contact pairs.

        Returns
        -------
        col1 : AtomGroup
            AtomGroup object for the first column of the contact pairs.
        """
        return self.atoms[:, 0]

    @property
    def col2(self):
        """Return the second column of the contact pairs.

        Returns
        -------
        col2 : AtomGroup
            AtomGroup object for the second column of the contact pairs.
        """
        return self.atoms[:, 1]

    def contacts(
        self,
        output: Literal["dict", "dataframe"],
        dftype: Literal["compact", "expanded", "print"] = "compact",
    ) -> Union[Dict[Any, Any], pd.DataFrame]:
        """Get the covalent contacts.

        Parameters
        ----------
        output : Literal["dict", "dataframe"], optional
            The output of the output, by default "dict"
        """
        if self.pairs is None or self.distances is None:
            raise ValueError(
                "No contacts found! Did you call the `compute_contacts` first?"
            )

        df = DataFrameFactory(self.atoms, self.distances, dftype).dataframe()
        if output == "dict":
            return df.to_dict(orient="list")
        elif output == "dataframe":
            return df
        else:
            raise ValueError(f"Invalid output type: {output}")

    def filter_distances(self):
        """Return distances between contact pairs."""
        indices = matching_indices(self.neighbors.pairs, self.pairs)
        return self.neighbors.distances[indices]

    @abstractmethod
    def compute_contacts(self, **kwargs) -> np.ndarray:
        """Main method for computing the contact pairs. Must be implemented by subclasses."""
