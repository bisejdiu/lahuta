"""The main entry point for users to access the various functionalities
of the `contacts` package. It imports all necessary functions for computing interactions, thus
providing an easy and concise way to use the package's functionalities.

Each function takes a `NeighborPairs` object, which represents precomputed neighbor relationships between atoms, 
and an optional distance cutoff parameter for certain types of contacts. Each function returns a new `NeighborPairs` 
object that represents the specific contacts computed by that function.

Notes:
    Each atomic contact type is computed based on a set of conditions, primarily involving the types of the two 
    interacting atoms and the distance between them. For example, carbonyl contacts are calculated by filtering 
    for neighbor pairs where one atom is a carbonyl oxygen and the other is a carbonyl carbon, and the distance 
    between them is less than a predefined cutoff.
    
    Most contact calculation functions also accept a `distance` argument specifying the cutoff distance for 
    considering two atoms as neighbors in the contact type calculation. These cutoffs are defined in the 
    `lahuta.config.defaults` module and can be adjusted as per user requirements.

Example:
    ``` py
    from lahuta.contacts import F

    # Define your system and compute neighbors
    luni = Luni(...)
    ns = luni.compute_neighbors()

    # Compute aromatic neighbors
    aromatic_nbs = F.aromatic_neighbors(ns)

    # Compute hydrophobic neighbors
    hydrophobic_nbs = F.hydrophobic_neighbors(ns)

    # Compute ionic neighbors
    ionic_nbs = F.ionic_neighbors(ns)

    # Compute plane-plane neighbors
    pp_nbs = F.plane_plane_neighbors(ns)
    ```

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
