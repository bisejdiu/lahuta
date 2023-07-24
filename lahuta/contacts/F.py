"""
Module: F.py

This module is designed as the main entry point for users to access the various functionalities
of the `contacts` package. It imports all necessary functions for computing interactions, thus
providing an easy and concise way to use the package's functionalities.

Functions:
    All functions for computing interactions between atoms and planes, and plane-plane interactions
    are available via this module.

Import Example:
    To use this module, you would typically import it as follows:

    import lahuta.contacts as F

    And then call the desired functions as methods of the imported module. For example, to compute 
    aromatic neighbors, you would do:

    F.aromatic_neighbors(ns)

    where `ns` is a NeighborPairs object.
    
Please refer to the individual function documentation for details about their usage and parameters.

Example:
    import lahuta.contacts as F

    # Define your system and compute neighbors
    universe = Universe(...)
    ns = universe.compute_neighbors()

    # Compute aromatic neighbors
    aromatic_nbs = F.aromatic_neighbors(ns)

    # Compute hydrophobic neighbors
    hydrophobic_nbs = F.hydrophobic_neighbors(ns)

    # Compute ionic neighbors
    ionic_nbs = F.ionic_neighbors(ns)

    # Compute plane-plane neighbors
    pp_nbs = F.plane_plane_neighbors(ns)

"""

from .atom_plane import carbon_pi, cation_pi, donor_pi, sulphur_pi
from .contacts import (
    aromatic_neighbors,
    carbonyl_neighbors,
    covalent_neighbors,
    hbond_neighbors,
    hydrophobic_neighbors,
    ionic_neighbors,
    metalic_neighbors,
    polar_hbond_neighbors,
    vdw_neighbors,
    weak_hbond_neighbors,
    weak_polar_hbond_neighbors,
)
from .plane_plane import plane_plane_neighbors

__all__ = [
    "aromatic_neighbors",
    "carbonyl_neighbors",
    "covalent_neighbors",
    "hbond_neighbors",
    "hydrophobic_neighbors",
    "ionic_neighbors",
    "metalic_neighbors",
    "polar_hbond_neighbors",
    "vdw_neighbors",
    "weak_hbond_neighbors",
    "weak_polar_hbond_neighbors",
    "donor_pi",
    "sulphur_pi",
    "carbon_pi",
    "cation_pi",
    "plane_plane_neighbors",
]
