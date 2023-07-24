"""
Module: carbon_pi.py

This module provides class-level API for computing carbon-pi contacts, based on atom-plane interactions. 

Warning: 
    Due to re-computation, using the `CarbonPi` class may result in slower performance. 
    We recommend using the `AtomPlaneContacts` class from the `atom_plane` module for 
    improved efficiency. 

Classes:
    CarbonPi: A class to compute carbon-pi contacts. When instantiated, a warning about 
    potential performance issues is emitted.

Usage:
    universe = Universe(...)
    ns = universe.compute_neighbors()
    carbon_pi = CarbonPi(ns)
    results = carbon_pi.compute().results
"""

from warnings import warn

from lahuta.core.neighbors import NeighborPairs

from .atom_plane import DEFAULT_CONTACT_DISTS, carbon_pi
from .base import ContactAnalysis


class CarbonPi(ContactAnalysis):
    distance = DEFAULT_CONTACT_DISTS["carbon_pi"]
    cache = False

    def __init__(self, ns: NeighborPairs):
        super().__init__(ns)
        warn(
            "Using the `CarbonPi` class may result in slower performance due to "
            "re-computation. Consider using the `AtomPlaneContacts` class from the "
            "`atom_plane` module for improved efficiency.",
            RuntimeWarning,
        )

    def compute(self) -> NeighborPairs:
        """Compute aromatic contacts."""
        return carbon_pi(self.ns, cache=self.cache)
