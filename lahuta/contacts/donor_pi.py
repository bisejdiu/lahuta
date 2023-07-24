"""
Module: donor_pi.py

This module provides class-level API for computing donor-pi contacts, based on atom-plane interactions. 

Warning: 
    Due to re-computation, using the `DonorPi` class may result in slower performance. 
    We recommend using the `AtomPlaneContacts` class from the `atom_plane` module for 
    improved efficiency. 

Classes:
    DonorPi: A class to compute donor-pi contacts. When instantiated, a warning about 
    potential performance issues is emitted.

Usage:
    universe = Universe(...)
    ns = universe.compute_neighbors()
    donor_pi = DonorPi(ns)
    results = donor_pi.compute().results
"""
from warnings import warn

from lahuta.core.neighbors import NeighborPairs

from .atom_plane import DEFAULT_CONTACT_DISTS, donor_pi
from .base import ContactAnalysis


class DonorPi(ContactAnalysis):
    distance = DEFAULT_CONTACT_DISTS["donor_pi"]
    cache = False

    def __init__(self, ns: NeighborPairs):
        super().__init__(ns)
        warn(
            "Using the `DonorPi` class may result in slower performance due to "
            "re-computation. Consider using the `AtomPlaneContacts` class from the "
            "`atom_plane` module for improved efficiency.",
            RuntimeWarning,
        )

    def compute(self) -> NeighborPairs:
        """Compute aromatic contacts."""
        return donor_pi(self.ns, cache=self.cache)
