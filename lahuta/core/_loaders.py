from abc import ABC, abstractmethod
from typing import Literal

import gemmi
import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.core.groups import AtomGroup

from lahuta.core.arc import ARC, Atoms, Chains, Residues
from lahuta.core.obmol import OBMol


class BaseLoader(ABC):
    def __init__(self, file_path):
        self.file_path = file_path
        self._chains = None
        self._residues = None
        self._atoms = None
        # self._coords_array = None

        self.structure = None
        self.ag = None
        self.arc = None

    # def _validate_access(self, attr_name):
    #     if getattr(self, attr_name) is None:
    #         raise ValueError(f"No {attr_name} in the loader")

    @property
    def n_atoms(self):
        # self._validate_access("_atoms")
        return len(self.atoms)  # type: ignore

    @property
    def chains(self):
        # self._validate_access("_chains")
        return self.arc.chains

    @property
    def residues(self):
        # self._validate_access("_residues")
        return self.arc.residues

    @property
    def atoms(self):
        # self._validate_access("_atoms")
        return self.arc.atoms

    # @property
    # def coords_array(self):
    #     return self._coords_array

    def to(self, object_type: Literal["mol", "mda"], *args, **kwargs):
        method_str = f"to_{object_type}"
        if hasattr(self, method_str):
            return getattr(self, method_str)()

        raise ValueError(f"Object type {object_type} is not supported")

    @abstractmethod
    def to_mda(self):
        ...

    @abstractmethod
    def to_mol(self):
        ...

    def __iter__(self):
        yield self.atoms
        yield self.residues
        yield self.chains


class GemmiLoader(BaseLoader):
    def __init__(self, file_path, is_pdb=False):
        super().__init__(file_path)
        if is_pdb:
            structure = gemmi.read_pdb(self.file_path)
            block = structure.make_mmcif_document().sole_block()
        else:
            block = gemmi.cif.read(self.file_path).sole_block()
            structure = gemmi.make_structure_from_block(block)

        self.structure = structure
        atom_site_data = block.get_mmcif_category("_atom_site.")

        self.arc = ARC(self, atom_site_data)
        # self._coords_array = self.extract_positions(atom_site_data)
        self.arc.atoms.coordinates = self.extract_positions(atom_site_data)

        self.ag = None

    def extract_positions(self, atom_site_data):
        coords_array = np.zeros((self.n_atoms, 3))
        coords_array[:, 0] = atom_site_data.get("Cartn_x")
        coords_array[:, 1] = atom_site_data.get("Cartn_y")
        coords_array[:, 2] = atom_site_data.get("Cartn_z")

        return coords_array

    def to_mda(self):
        if self.ag is not None:
            return self.ag

        # Create a structured array to ensure unique values for each combination of resname, resid, and chain_id
        struct_arr = np.rec.fromarrays(
            [self.arc.residues.resnames, self.arc.residues.resids, self.arc.chains.ids],
            names="resnames, resids, chain_ids",
        )

        # Use factorize to get the labels and unique values
        resindices, uniques = pd.factorize(struct_arr)

        resnames, resids, chain_ids = (
            uniques["resnames"],
            uniques["resids"],
            uniques["chain_ids"],
        )

        # Create a new Universe
        uv = mda.Universe.empty(
            n_atoms=self.arc.atoms.ids.size,
            n_residues=uniques.size,
            n_segments=chain_ids.size,
            atom_resindex=resindices,
            residue_segindex=chain_ids,
            trajectory=True,
        )

        # Add topology attributes
        uv.add_TopologyAttr("names", self.arc.atoms.names)
        uv.add_TopologyAttr("type", self.arc.atoms.types)
        uv.add_TopologyAttr("elements", self.arc.atoms.elements)
        uv.add_TopologyAttr("resnames", resnames)
        uv.add_TopologyAttr("resids", resids)
        uv.add_TopologyAttr("segids", chain_ids)

        uv.atoms.positions = self.arc.atoms.coordinates  # type: ignore

        self.ag = uv.atoms

        return self.ag

    def to_mol(self):
        obmol = OBMol()
        obmol.create_mol(
            self.arc,
            # self.coords_array,
            self.structure.connections,
        )

        return obmol.mol


class TopologyLoader(BaseLoader):
    def __init__(self, *paths):
        print("Using TopologyLoader")
        file_path = paths[0]
        super().__init__(file_path)
        universe = mda.Universe(self.file_path)
        self.ag: AtomGroup = universe.atoms  # type: ignore
        assert self.ag is not None
        if len(paths) > 1:
            self.ag.universe.load_new(paths[1:], format=None, in_memory=False)

        self.arc = ARC(self, self.ag)  # positions are set when using mda.Universe

    def to_mda(self):
        return self.ag

    def to_mol(self):
        obmol = OBMol()
        obmol.create_mol(
            self.arc,
            # self.coords_array,
            self.structure,
        )

        return obmol.mol

    @classmethod
    def from_mda(cls, ag):
        top_loader = cls.__new__(cls)
        top_loader.ag = ag.copy()
        top_loader.ag._u = ag.universe.copy()
        top_loader.structure = None

        top_loader.arc = ARC(top_loader, top_loader.ag)

        return top_loader
