"""
Placeholder for the universe module.
"""

import lahuta.contacts as F
from lahuta.config.defaults import CONTACTS

from .base import ContactAnalysis


class HydrophobicContacts(ContactAnalysis):
    """A class to compute hydrophobic contacts."""

    distance = CONTACTS["hydrophobic"]["distance"]

    def compute(self):
        """Compute hydrophobic contacts."""
        return F.hydrophobic_neighbors(self.ns, self.distance)
