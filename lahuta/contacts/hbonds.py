"""
Module: hbonds.py

This module defines a class for computing hbond-type contacts using a class-based approach. 
The HBondContacts class inherits from the base ContactAnalysis class and 
implements the `compute` method for hbond contact computation.

Class:
    HBondContacts(ContactAnalysis): Computes hbond contacts.
    WeakHBondContacts(ContactAnalysis): Computes weak hbond contacts.
    PolarHBondContacts(ContactAnalysis): Computes polar hbond contacts.
    WeakPolarHBondContacts(ContactAnalysis): Computes weak polar hbond contacts.
                                       
Example:
    universe = Universe(...)
    ns = universe.compute_neighbors()

    hbonds = HBondContacts(ns)
    print(hbonds.results)
"""

import lahuta.contacts as F
from lahuta.config.defaults import CONTACTS
from lahuta.core.neighbors import NeighborPairs

from .base import ContactAnalysis


class HBondContacts(ContactAnalysis):
    """A class to find hbond contacts between atoms in a molecule."""

    def compute(self) -> NeighborPairs:
        """Compute hbond contacts."""
        return F.hbond_neighbors(self.ns)


class WeakHBondContacts(ContactAnalysis):
    """A class to find weak hbond contacts between atoms in a molecule."""

    def compute(self) -> NeighborPairs:
        """Compute weak hbond contacts."""
        return F.weak_hbond_neighbors(self.ns)


class PolarHBondContacts(ContactAnalysis):
    """A class to find polar hbond contacts between atoms in a molecule."""

    distance = CONTACTS["hbond"]["polar distance"]

    def compute(self) -> NeighborPairs:
        """Compute polar hbond contacts."""
        return F.polar_hbond_neighbors(self.ns, self.distance)


class WeakPolarHBondContacts(ContactAnalysis):
    """A class to find weak polar hbond contacts between atoms in a molecule."""

    distance = CONTACTS["weak hbond"]["weak polar distance"]

    def compute(self) -> NeighborPairs:
        """Compute weak polar hbond contacts."""
        return F.weak_polar_hbond_neighbors(self.ns, self.distance)
