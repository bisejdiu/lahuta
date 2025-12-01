"""
Backbone for a inter-operability layer between Lahuta and other libraries.
Currently, this is more an ambition than a reality.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..lib import lahuta as _lib


# fmt: off
@dataclass
class MolecularSystem:
    """High-level representation of a molecular system, useful for interfacing with other libraries."""

    atom_indices:   list[int]
    atomic_numbers: list[int]
    atom_names:     list[str]
    resids:         list[int]
    resnames:       list[str]
    chainlabels:    list[str]
    positions:      np.ndarray
    path:           Path | None

    @classmethod
    def from_file(cls, path: str | Path) -> "MolecularSystem":
        """Create MolecularSystem from structure file."""
        temp_luni = _lib.LahutaSystem(str(path))

        return cls(
            atom_indices   = temp_luni.props.indices.astype(int).tolist(),
            atomic_numbers = temp_luni.props.atom_nums.astype(int).tolist(),
            atom_names     = temp_luni.props.names.tolist(),
            resids         = temp_luni.props.resids.astype(int).tolist(),
            resnames       = temp_luni.props.resnames.tolist(),
            chainlabels    = temp_luni.props.chainlabels.tolist(),
            positions      = temp_luni.props.positions,
            path=Path(path) if isinstance(path, str) else path,
        )

    def from_subset(self, indices: list[int]) -> "MolecularSystem":
        """Create new system containing only specified atoms."""
        mask = np.array([i in indices for i in range(len(self.atom_indices))])

        return MolecularSystem(
            atom_indices    = [self.atom_indices[i]   for i in indices],
            atomic_numbers  = [self.atomic_numbers[i] for i in indices],
            atom_names      = [self.atom_names[i]     for i in indices],
            resids          = [self.resids[i]         for i in indices],
            resnames        = [self.resnames[i]       for i in indices],
            chainlabels     = [self.chainlabels[i]    for i in indices],
            positions=self.positions[mask],
            path=self.path,
        )

    @classmethod
    def from_ir(cls, ir: _lib.IR) -> "MolecularSystem":
        """Create MolecularSystem from intermediate representation."""
        return cls(
            atom_indices   = ir.atom_indices,
            atomic_numbers = ir.atomic_numbers,
            atom_names     = ir.atom_names,
            resids         = ir.resids,
            resnames       = ir.resnames,
            chainlabels    = ir.chainlabels,
            positions      = np.array(ir.positions),
            path=None,
        )

    # NOTE: some of these may be handled directly by lhxx
    def to_ir(self) -> _lib.IR:
        """Convert to intermediate representation for C++ bindings."""
        return _lib.IR(
            self.atom_indices,
            self.atomic_numbers,
            self.atom_names,
            self.resids,
            self.resnames,
            self.chainlabels,
            self.positions.tolist(),
        )

    def to_mda(self):
        """Convert to MDAnalysis Universe."""
        pass

    def to_rdkit(self):
        """Convert to RDKit Mol object."""
        pass

    def to_openmm(self):
        """Convert to OpenMM Topology and Positions."""
        pass

    def to_openbabel(self):
        """Convert to OpenBabel OBMol."""
        pass
