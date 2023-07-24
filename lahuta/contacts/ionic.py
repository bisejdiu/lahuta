"""
Module: ionic.py

This module defines a class for computing ionic contacts using a class-based approach. 
The IonicContacts class inherits from the base ContactAnalysis class and 
implements the `compute` method for ionic contact computation.

Class:
    IonicContacts(ContactAnalysis): Computes ionic contacts.
                                       
Example:
    universe = Universe(...)
    ns = universe.compute_neighbors()

    ionic = IonicContacts(ns)
    print(ionic.results)
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
