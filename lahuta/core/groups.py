"""
Placeholder for the groups module.
"""

import MDAnalysis as mda
import numpy as np
from MDAnalysis.core.topologyattrs import AtomAttr
from MDAnalysis.lib.nsgrid import FastNS

# import from local utils module
from ..utils.mda import mda_psuedobox_from_atomgroup
from .neighbors import NeighborPairs


class AtomGroup(mda.AtomGroup):
    """A subclass of the MDAnalysis AtomGroup class."""

    def __init__(self, *args, **kwargs):
        """Initialize the AtomGroup."""
        super().__init__(*args, **kwargs)

    def compute_neighbors(
        self,
        radius=5.0,
        ignore_hydrogens=True,
        # hbonds=False,
        skip_adjacent=True,
        res_dif=1,
    ):
        """Compute the neighbors of each atom in the Universe.

        Parameters
        ----------
        radius : float, optional
            The cutoff radius. Default is 5.0.
        ignore_hydrogens : bool, optional
            Whether to ignore hydrogens. Default is True.
        hbonds : bool, optional
            Whether to include hydrogen bonds. Default is False.

        Returns
        -------
        neighbors : np.ndarray
            An array of shape (n_atoms, n_neighbors) where each row contains the
            indices of the neighbors of the atom in the row.
        """
        if ignore_hydrogens:
            atomgroup = self.select_atoms("not name H*")
        else:
            atomgroup = self.atoms

        pairs, distances = self.get_neighbors(atomgroup, radius)

        if skip_adjacent:
            idx = self._remove_adjacent_residue_pairs(pairs, res_dif=res_dif)
            pairs = pairs[idx]
            distances = distances[idx]

        return NeighborPairs(self.atoms, pairs, distances)

    def get_neighbors(self, atomgroup=None, radius=5.0):
        """Get the neighbors of an atomgroup.

        Parameters
        ----------
        atomgroup : AtomGroup
            The atomgroup to get the neighbors of.

        Returns
        -------
        pairs : np.ndarray
            An array of shape (n_pairs, 2) where each row contains the indices of
            the atoms in the pair.
        distances : np.ndarray
            An array of shape (n_pairs,) containing the distances of each pair.
        """
        if atomgroup is None:
            atomgroup = self.atoms

        shift_coords, pseudobox = mda_psuedobox_from_atomgroup(atomgroup)

        gridsearch = FastNS(
            cutoff=radius, coords=shift_coords, box=pseudobox, pbc=False
        )
        neighbors = gridsearch.self_search()

        return atomgroup[neighbors.get_pairs()].indices, neighbors.get_pair_distances()

    def _remove_adjacent_residue_pairs(self, pairs, res_dif=1):
        """Remove pairs where the difference in residue ids is less than `res_dif`.

        Parameters
        ----------
        uniag : AtomGroup/Universe
            The atomgroup containing the atoms in the pairs (can also be a Universe).
        pairs : np.ndarray
            An array of shape (n_pairs, 2) where each row is a pair of atom indices.
        res_dif : int, optional
            The difference in residue ids to remove. Default is 1.

        Returns
        -------
        pairs : np.ndarray
            An array of shape (n_pairs, 2) where each row is a pair of atom indices.

        """
        # so it works for both AtomGroup and Universe
        resids = self.atoms.resids[pairs]
        return np.any(np.abs(resids - resids[:, ::-1]) > res_dif, axis=1)

    def __str__(self) -> str:
        return f"<Lahuta AtomGroup containing {self.indices.size} atoms>"

    def __repr__(self) -> str:
        return self.__str__()


class VdWRadii(AtomAttr):
    """VdW radii for all atoms in the Universe."""

    attrname = "vdw_radii"
    singular = "vdw_radius"

    @staticmethod
    def _gen_initial_values(n_atoms, n_residues, n_segments, radii):
        """Generate the initial values for the VdW radii."""
        return radii


class CovRadii(AtomAttr):
    """Covalent radii for all atoms in the Universe."""

    attrname = "cov_radii"
    singular = "cov_radius"

    @staticmethod
    def _gen_initial_values(n_atoms, n_residues, n_segments, radii):
        """Generate the initial values for the VdW radii."""
        return radii[:, 1]


class AtomTypes(AtomAttr):
    """Atom types for all atoms in the Universe."""

    attrname = "atom_types"
    singular = "atom_type"

    @staticmethod
    def _gen_initial_values(n_atoms, n_residues, n_segments, atypes):
        """Generate the initial values for the atom types."""
        return atypes
