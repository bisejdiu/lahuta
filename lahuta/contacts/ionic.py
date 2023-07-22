"""
Placeholder for the universe module.
"""

import lahuta.contacts as F
from lahuta.config.defaults import CONTACTS
from lahuta.core.neighbors import NeighborPairs

from .base import ContactAnalysis


class IonicContacts(ContactAnalysis):
    """A class to find ionic contacts between atoms in a molecule."""

    distance = CONTACTS["ionic"]["distance"]

    def compute(self) -> NeighborPairs:
        """Compute ionic contacts."""
        return F.ionic_neighbors(self.ns, self.distance)
