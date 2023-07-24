"""
Module: covalent.py

This module defines a class for computing covalent contacts using a class-based approach. 
The CovalentContacts class inherits from the base ContactAnalysis class and 
implements the `compute` method for covalent contact computation.

Class:
    CovalentContacts(ContactAnalysis): Computes covalent contacts.


Note:
    This contact is slightly slower because it requires a loop over all atoms in the molecule 
    to check for pairs that are covalently bonded. 
                               
Example:
    universe = Universe(...)
    ns = universe.compute_neighbors()

    cov = CovalentContacts(ns)
    print(cov.results)
"""


import lahuta.contacts as F
from lahuta.core.neighbors import NeighborPairs

from .base import ContactAnalysis


class CovalentContacts(ContactAnalysis):
    """A class to compute covalent contacts between atoms in a molecule."""

    def compute(self) -> NeighborPairs:
        """Compute covalent contacts."""
        return F.covalent_neighbors(self.ns)
