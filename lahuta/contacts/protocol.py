"""
Placeholder
"""

from typing_extensions import Protocol

import numpy as np
import MDAnalysis as mda

from ..core.neighbors import NeighborPairs


class ContactStrategy(Protocol):
    """A protocol for contact strategies."""

    universe: mda.Universe
    neighbor_pairs: NeighborPairs

    @property
    def contact_pairs(self, **kwargs) -> np.ndarray:
        """Retrieve the contact pairs indices."""

    @property
    def contact_atoms(self, **kwargs) -> mda.AtomGroup:
        """Retrieve the contact atoms."""

    def _contact_pairs(self, **kwargs) -> np.ndarray:
        """Main method for computing the contact pairs."""
