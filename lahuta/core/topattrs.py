import numpy as np
from MDAnalysis.core.topologyattrs import AtomAttr


class VdWRadii(AtomAttr):
    """VdW radii for all atoms in the Universe."""

    attrname = "vdw_radii"
    singular = "vdw_radius"

    # @staticmethod
    # def _gen_initial_values(n_atoms, n_residues, n_segments):
    #     """Generate the initial values for the VdW radii."""
    #     return np.zeros(n_atoms)


class CovRadii(AtomAttr):
    """Covalent radii for all atoms in the Universe."""

    attrname = "cov_radii"
    singular = "cov_radius"

    # @staticmethod
    # def _gen_initial_values(n_atoms, n_residues, n_segments, radii):
    #     """Generate the initial values for the VdW radii."""
    #     return radii[:, 1]


class AtomTypes(AtomAttr):
    """Atom types for all atoms in the Universe."""

    attrname = "atom_types"
    singular = "atom_type"

    # @staticmethod
    # def _gen_initial_values(n_atoms, n_residues, n_segments):
    #     """Generate the initial values for the atom types."""
    #     return np.zeros(n_atoms, dtype=np.int32)
