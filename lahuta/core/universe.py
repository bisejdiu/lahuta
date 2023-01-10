"""
Placeholder for the universe module.
"""

from pathlib import Path

import MDAnalysis as mda
import numpy as np

from ..utils.atom_types import (
    assign_atom_types,
    assign_radii,
    find_hydrogen_bonded_atoms,
)
from .groups import AtomGroup
from .obabel import OBMol


class Universe(mda.Universe):
    """A subclass of MDAnalysis Universe that adds additional attributes and methods."""

    def __init__(self, *args, **kwargs):
        """A subclass of MDAnalysis Universe that adds some extra functionality."""
        # get the first argument from the args list
        self.filename = args[0]

        suffix = "pdb"
        if isinstance(self.filename, str):
            suffix = self.filename.split(".")[-1]
        elif isinstance(self.filename, Path):
            suffix = self.filename.suffix[1:]
            self.filename = str(self.filename)

        if suffix == "cif":
            from io import StringIO

            from lahuta.utils.cif_converter import convert_cif_to_pdb, read_cif_as_pdb

            # pdb_str, mol = convert_cif_to_pdb(self.filename)
            pdb_str, mol = read_cif_as_pdb(self.filename)
            super().__init__(
                mda.lib.util.NamedStream(StringIO(pdb_str), "dummy.pdb"), **kwargs
            )
            self.mol = mol
        elif suffix == "pdb":

            super().__init__(*args, **kwargs)
            self.mol = OBMol(str(self.filename))

        self.atoms = AtomGroup(self.atoms)
        self.hbond_array = find_hydrogen_bonded_atoms(self.mol)

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

    def select_atoms(self, *args, **kwargs) -> AtomGroup:
        """Select atoms.

        Wrapper around MDAnalysis Universe.select_atoms that returns a LahutaAtomGroup.
        """
        return self.atoms.select_atoms(*args, **kwargs)  # type: ignore

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
