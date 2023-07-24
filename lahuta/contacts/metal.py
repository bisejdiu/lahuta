"""
Module: metal.py

This module defines a class for computing metal contacts using a class-based approach. 
The MetalicContacts class inherits from the base ContactAnalysis class and 
implements the `compute` method for metal contact computation.

Class:
    MetalicContacts(ContactAnalysis): Computes metal contacts.
                                       
Example:
    universe = Universe(...)
    ns = universe.compute_neighbors()

    metal_contacts = MetalicContacts(ns)
    print(metal_contacts.results)
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
