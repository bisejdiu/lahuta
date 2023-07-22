"""
Placeholder for the universe module.
"""

import lahuta.contacts as F
from lahuta.config.defaults import CONTACTS
from lahuta.core.neighbors import NeighborPairs

# from ..config.atoms import METALS
from .base import ContactAnalysis


class MetalicContacts(ContactAnalysis):
    """A class to compute metal contacts."""

    distance = CONTACTS["metal"]["distance"]

    def compute(self) -> NeighborPairs:
        """Compute metal contacts."""
        return F.metalic_neighbors(self.ns, self.distance)
