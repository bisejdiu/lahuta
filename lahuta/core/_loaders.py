from abc import ABC, abstractmethod
from typing import Literal

import gemmi
import MDAnalysis as mda
import numpy as np
from MDAnalysis.core.groups import AtomGroup

from lahuta.core.cra import Atoms, Chains, Residues
from lahuta.core.obmol import OBMol


class BaseLoader(ABC):
    def __init__(self, file_path):
        self.file_path = file_path
        self._chains = None
        self._residues = None
        self._atoms = None
        self._coords_array = None

        self.structure = None
        self.ag = None

    def _validate_access(self, attr_name):
        if getattr(self, attr_name) is None:
            raise ValueError(f"No {attr_name} in the loader")

    @property
    def n_atoms(self):
        self._validate_access("_atoms")
        return len(self._atoms)  # type: ignore

    @property
    def chains(self):
        self._validate_access("_chains")
        return self._chains

    @property
    def residues(self):
        self._validate_access("_residues")
        return self._residues

    @property
    def atoms(self):
        self._validate_access("_atoms")
        return self._atoms

    @property
    def coords_array(self):
        return self._coords_array

    def to(self, object_type: Literal["mol", "mda"], *args, **kwargs):
        method_str = f"to_{object_type}"
        if hasattr(self, method_str):
            return getattr(self, method_str)()

        raise ValueError(f"Object type {object_type} is not supported")

    @abstractmethod
    def create(self):
        ...

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
        # print("Using GemmiLoader")
        super().__init__(file_path)
        if is_pdb:
            structure = gemmi.read_pdb(self.file_path)
            block = structure.make_mmcif_document().sole_block()
        else:
            block = gemmi.cif.read(self.file_path).sole_block()
            structure = gemmi.make_structure_from_block(block)

        self.structure = structure
        atom_site_data = block.get_mmcif_category("_atom_site.")
        self._chains, self._residues, self._atoms = self.create(atom_site_data)
        self._coords_array = self.extract_positions(atom_site_data)

    def create(self, atom_site_data):
        chains = Chains().from_gemmi(atom_site_data)
        residues = Residues().from_gemmi(atom_site_data)
        atoms = Atoms().from_gemmi(atom_site_data)
        return chains, residues, atoms

    def extract_positions(self, atom_site_data):
        coords_array = np.zeros((self.n_atoms, 3))
        coords_array[:, 0] = atom_site_data.get("Cartn_x")
        coords_array[:, 1] = atom_site_data.get("Cartn_y")
        coords_array[:, 2] = atom_site_data.get("Cartn_z")

        return coords_array

    def to_mda(self):
        if self.ag is not None:
            return self.ag

        # TODO: Add icodes and ids to the universe
        resnames, resids, chain_ids = [], [], []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    resids.append(residue.seqid.num)
                    resnames.append(residue.name)
                    chain_ids.append(self.chains.mapping[chain.name])

        universe = mda.Universe.empty(
            n_atoms=self.n_atoms,
            n_residues=len(resids),
            atom_resindex=self.residues.resindices,
            residue_segindex=chain_ids,
            trajectory=True,
        )

        universe.add_TopologyAttr("names", self.atoms.names)
        universe.add_TopologyAttr("type", self.atoms.types)
        universe.add_TopologyAttr("elements", self.atoms.elements)
        universe.add_TopologyAttr("resnames", resnames)
        universe.add_TopologyAttr("resids", resids)
        universe.add_TopologyAttr("segids", np.array(["PROT"], dtype=object))

        universe.atoms.positions = self.coords_array  # type: ignore

        self.ag = universe.atoms

        return self.ag

    def to_mol(self):
        obmol = OBMol()
        obmol.create_mol(
            (self.chains, self.residues, self.atoms),
            self.coords_array,
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
        self._chains, self._residues, self._atoms = self.create()
        self._coords_array = self.ag.atoms.positions  # type: ignore

    def create(self):
        chains = Chains().from_mda(self.ag)
        residues = Residues().from_mda(self.ag)
        atoms = Atoms().from_mda(self.ag)
        return chains, residues, atoms

    def to_mda(self):
        return self.ag

    def to_mol(self):
        obmol = OBMol()
        obmol.create_mol(
            (self.chains, self.residues, self.atoms),
            self.coords_array,
            self.structure,
        )

        return obmol.mol

    @classmethod
    def from_mda(cls, ag):
        top_loader = cls.__new__(cls)
        top_loader.ag = ag.copy()
        top_loader.ag._u = ag.universe.copy()
        (
            top_loader._chains,
            top_loader._residues,
            top_loader._atoms,
        ) = top_loader.create()
        top_loader._coords_array = top_loader.ag.positions
        top_loader.structure = None

        return top_loader
