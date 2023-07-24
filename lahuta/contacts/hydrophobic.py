"""
Module: hydrophobic.py

This module defines a class for computing hydrophobic contacts using a class-based approach. 
The HydrophobicContacts class inherits from the base ContactAnalysis class and 
implements the `compute` method for hydrophobic contact computation.

Class:
    HydrophobicContacts(ContactAnalysis): Computes hydrophobic contacts.
                                       
Example:
    universe = Universe(...)
    ns = universe.compute_neighbors()

    hbc = HydrophobicContacts(ns)
    print(hbc.results)
"""

import lahuta.contacts as F
from lahuta.config.defaults import CONTACTS
from lahuta.core.neighbors import NeighborPairs

from .base import ContactAnalysis


class HydrophobicContacts(ContactAnalysis):
    """A class to compute hydrophobic contacts."""

    distance = CONTACTS["hydrophobic"]["distance"]

    def compute(self) -> NeighborPairs:
        """Compute hydrophobic contacts."""
        return F.hydrophobic_neighbors(self.ns, self.distance)
