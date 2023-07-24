"""
Module: cation_pi.py

This module provides class-level API for computing cation-pi contacts, based on atom-plane interactions. 

Warning: 
    Due to re-computation, using the `CationPi` class may result in slower performance. 
    We recommend using the `AtomPlaneContacts` class from the `atom_plane` module for 
    improved efficiency. 

Classes:
    CationPi: A class to compute cation-pi contacts. When instantiated, a warning about 
    potential performance issues is emitted.

Usage:
    universe = Universe(...)
    ns = universe.compute_neighbors()
    cation_pi = CationPi(ns)
    results = cation_pi.compute().results
"""
from warnings import warn

from lahuta.core.neighbors import NeighborPairs

from .atom_plane import DEFAULT_CONTACT_DISTS, cation_pi
from .base import ContactAnalysis


class CationPi(ContactAnalysis):
    distance = DEFAULT_CONTACT_DISTS["cation_pi"]
    cache = False

    def __init__(self, ns: NeighborPairs):
        super().__init__(ns)
        warn(
            "Using the `CationPi` class may result in slower performance due to "
            "re-computation. Consider using the `AtomPlaneContacts` class from the "
            "`atom_plane` module for improved efficiency.",
            RuntimeWarning,
        )

    def compute(self) -> NeighborPairs:
        """Compute aromatic contacts."""
        return cation_pi(self.ns, cache=self.cache)
