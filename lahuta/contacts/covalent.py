"""
Placeholder for the universe module.
"""

import lahuta.contacts as F
from lahuta.core.neighbors import NeighborPairs

from .base import ContactAnalysis


class CovalentContacts(ContactAnalysis):
    """A class to compute covalent contacts between atoms in a molecule."""

    def compute(self) -> NeighborPairs:
        """Compute covalent contacts."""
        return F.covalent_neighbors(self.ns)
