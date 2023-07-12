import MDAnalysis as mda
import numpy as np

from lahuta.core.neighbor_finder import NeighborSearch


class AtomGroup(mda.AtomGroup):
    """A subclass of the MDAnalysis AtomGroup class."""

    def __init__(self, *args, **kwargs):
        """Initialize the AtomGroup."""
        super().__init__(*args, **kwargs)

    def compute_neighbors(
        self,
        radius=5.0,
        ignore_hydrogens=True,
        skip_adjacent=True,
        res_dif=1,
    ):
        """
        Compute the neighbors of each atom in the Universe.

        Args:
        ----
        radius (float, optional): The cutoff radius. Default is 5.0.
        ignore_hydrogens (bool, optional): Whether to ignore hydrogens. Default is True.
        skip_adjacent (bool, optional): Whether to skip adjacent. Default is True.
        res_dif (int, optional): The residue difference to consider. Default is 1.

        Returns:
        -------
        NeighborPairs : np.ndarray
            An array of shape (n_atoms, n_neighbors) where each row contains the indices of the neighbors of the atom in the row.

        """
        if not self._ready:
            self.ready()

        neighbors = NeighborSearch(self)
        return neighbors.compute(
            radius=radius,
            ignore_hydrogens=ignore_hydrogens,
            skip_adjacent=skip_adjacent,
            res_dif=res_dif,
        )

    def __str__(self) -> str:
        return f"<Lahuta AtomGroup containing {self.indices.size} atoms>"

    def __repr__(self) -> str:
        return self.__str__()
