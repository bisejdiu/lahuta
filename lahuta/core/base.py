"""Base class for finding neighbors."""
from typing import Any

import numpy as np
from numpy.typing import NDArray

from lahuta._types.mdanalysis import AtomGroupType

IndexPairs = NDArray[np.int32]
Distances = NDArray[np.float32]
PairsDistances = tuple[IndexPairs, Distances]


class BaseNeighborSearch:
    """Base class for finding neighbors.

    The class provides methods to find neighbors of each atom in the universe and to remove pairs of atoms
    that are adjacent in the sequence, among other things.

    Args:
        mda (AtomGroupType): The AtomGroup containing the atoms.

    Attributes:
        ag_no_h (AtomGroup): Atom group of a universe excluding hydrogen atoms.
        og_resids (np.ndarray): The residue IDs of each atom in the universe.
    """

    def __init__(self, mda: AtomGroupType) -> None:
        self.ag_no_h = mda.select_atoms("not name H*")
        self.og_resids = mda.universe.atoms.resids
        self.chain_ids = mda.universe.atoms.chainIDs

    def compute(self, *args: Any, **kwargs: Any) -> PairsDistances: # noqa: ANN401
        """Compute the neighbors of each atom in the Universe."""
        raise NotImplementedError("This method should be overridden in subclasses")