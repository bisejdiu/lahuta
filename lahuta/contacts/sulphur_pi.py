"""
Module: sulphur_pi.py

This module provides class-level API for computing sulphur-pi contacts, based on atom-plane interactions. 

Warning: 
    Due to re-computation, using the `SulphurPi` class may result in slower performance. 
    We recommend using the `AtomPlaneContacts` class from the `atom_plane` module for 
    improved efficiency. 

Classes:
    SulphurPi: A class to compute sulphur-pi contacts. When instantiated, a warning about 
    potential performance issues is emitted.

Usage:
    universe = Universe(...)
    ns = universe.compute_neighbors()
    sulphur_pi = SulphurPi(ns)
    results = sulphur_pi.compute().results
"""
from warnings import warn

from lahuta.core.neighbors import NeighborPairs

from .atom_plane import DEFAULT_CONTACT_DISTS, sulphur_pi
from .base import ContactAnalysis


class SulphurPi(ContactAnalysis):
    distance = DEFAULT_CONTACT_DISTS["sulphur_pi"]
    cache = False

    def __init__(self, ns: NeighborPairs):
        super().__init__(ns)
        warn(
            "Using the `SulphurPi` class may result in slower performance due to "
            "re-computation. Consider using the `AtomPlaneContacts` class from the "
            "`atom_plane` module for improved efficiency.",
            RuntimeWarning,
        )

    def compute(self) -> NeighborPairs:
        """Compute aromatic contacts."""
        return sulphur_pi(self.ns, cache=self.cache)
