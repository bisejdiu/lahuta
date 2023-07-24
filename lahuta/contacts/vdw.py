"""
Module: vdw.py

This module defines a class for computing vanderwaals contacts using a class-based approach. 
The VanDerWaalsContacts class inherits from the base ContactAnalysis class and 
implements the `compute` method for vanderwaals contact computation.

Class:
    VanDerWaalsContacts(ContactAnalysis): Computes vanderwaals contacts.
                                       
Example:
    universe = Universe(...)
    ns = universe.compute_neighbors()

    vdw = VanDerWaalsContacts(ns)
    print(vdw.results)
"""

import lahuta.contacts as F
from lahuta.core.neighbors import NeighborPairs

from .base import ContactAnalysis


class VanDerWaalsContacts(ContactAnalysis):
    """A class to compute VdW contacts."""

    vdw_comp_factor = 0.1
    remove_clashes = True

    def compute(self) -> NeighborPairs:
        return F.vdw_neighbors(self.ns, self.vdw_comp_factor, self.remove_clashes)
