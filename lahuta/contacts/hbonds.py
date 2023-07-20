"""
Placeholder for the universe module.
"""

import lahuta.contacts as F
from lahuta.config.defaults import CONTACTS

from .base import ContactAnalysis


class HBondContacts(ContactAnalysis):
    """A class to find hbond contacts between atoms in a molecule."""

    def compute(self):
        """Compute hbond contacts."""
        return F.hbond_neighbors(self.ns)


class WeakHBondContacts(ContactAnalysis):
    """A class to find weak hbond contacts between atoms in a molecule."""

    def compute(self):
        """Compute weak hbond contacts."""
        return F.weak_hbond_neighbors(self.ns)


class PolarHBondContacts(ContactAnalysis):
    """A class to find polar hbond contacts between atoms in a molecule."""

    distance = CONTACTS["hbond"]["polar distance"]

    def compute(self):
        """Compute polar hbond contacts."""
        return F.polar_hbond_neighbors(self.ns, self.distance)


class WeakPolarHBondContacts(ContactAnalysis):
    """A class to find weak polar hbond contacts between atoms in a molecule."""

    distance = CONTACTS["weak hbond"]["weak polar distance"]

    def compute(self):
        """Compute weak polar hbond contacts."""
        return F.weak_polar_hbond_neighbors(self.ns, self.distance)
