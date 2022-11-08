"""
Placeholder for the universe module.
"""

import numpy as np
import MDAnalysis as mda

from .groups import AtomGroup

from .obabel import OBMol
from ..utils.atom_types import assign_atom_types, assign_radii


class Universe(mda.Universe):
    """A subclass of MDAnalysis Universe that adds additional attributes and methods."""

    def __init__(self, *args, **kwargs):
        """A subclass of MDAnalysis Universe that adds some extra functionality."""
        super().__init__(*args, **kwargs)

        self.mol = OBMol(self.filename)
        self.atoms = AtomGroup(self.atoms)

        self._extend_topology("vdw_radii", assign_radii(self.mol))
        self._extend_topology("atom_types", assign_atom_types(self.mol, self.atoms))

    def _extend_topology(self, attrname: str, values: np.ndarray):
        """Extend the topology with a new attribute.

        Parameters
        ----------
        attrname : str
            The name of the attribute.

        values : np.ndarray
            The values of the attribute.
        """
        self.add_TopologyAttr(attrname, values)

    def select_atoms(self, *args, **kwargs):
        """Select atoms.

        Wrapper around MDAnalysis Universe.select_atoms that returns a LahutaAtomGroup.
        """
        return self.atoms.select_atoms(*args, **kwargs)

    def compute_neighbors(self, *args, **kwargs):
        """Compute neighbors.
        See Also
        --------
        :meth:`Lahuta.core.groups.AtomGroup.compute_neighbors`
        """
        return self.atoms.compute_neighbors(*args, **kwargs)

    def __repr__(self):
        return f"<Lahuta Universe with {self.atoms.n_atoms} atoms>"

    def __str__(self):
        return self.__repr__()
