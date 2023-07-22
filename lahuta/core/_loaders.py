from abc import ABC, abstractmethod
from typing import Any, Dict, Iterator, Literal, Optional, Union, overload

import gemmi
import MDAnalysis as mda
import numpy as np
import pandas as pd
from numpy.typing import NDArray

from lahuta.core.arc import ARC, Atoms, Chains, Residues
from lahuta.core.obmol import OBMol
from lahuta.types.mdanalysis import AtomGroupType, UniverseType
from lahuta.types.openbabel import MolType


class BaseLoader(ABC):
    def __init__(self, file_path: str):
        self.file_path = file_path
        self._chains = None
        self._residues = None
        self._atoms = None

        self.structure = None
        # self.ag = None
        self.arc: Optional[ARC] = None

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def chains(self) -> Chains:
        if self.arc is None:
            raise ValueError("arc has not been initialized")
        return self.arc.chains

    @property
    def residues(self) -> Residues:
        if self.arc is None:
            raise ValueError("arc has not been initialized")
        return self.arc.residues

    @property
    def atoms(self) -> Atoms:
        if self.arc is None:
            raise ValueError("arc has not been initialized")
        return self.arc.atoms

    @overload
    def to(self, fmt: Literal["mda"]) -> AtomGroupType:
        ...

    @overload
    def to(self, fmt: Literal["mol"]) -> MolType:
        ...

    def to(self, fmt: Literal["mol", "mda"]) -> Union[MolType, AtomGroupType]:
        method_str = f"to_{fmt}"
        if hasattr(self, method_str):
            return getattr(self, method_str)()  # type: ignore

        raise ValueError(f"Object type {fmt} is not supported")

    @abstractmethod
    def to_mda(self) -> AtomGroupType:
        ...

    @abstractmethod
    def to_mol(self) -> MolType:
        ...

    def __iter__(self) -> Iterator[Union[Atoms, Residues, Chains]]:
        yield self.atoms
        yield self.residues
        yield self.chains


class GemmiLoader(BaseLoader):
    def __init__(self, file_path: str, is_pdb: bool = False):
        super().__init__(file_path)
        if is_pdb:
            structure: Any = gemmi.read_pdb(self.file_path)  # type: ignore
            block: Any = structure.make_mmcif_document().sole_block()
        else:
            block: Any = gemmi.cif.read(self.file_path).sole_block()  # type: ignore
            structure: Any = gemmi.make_structure_from_block(block)  # type: ignore

        self.structure = structure
        atom_site_data: Dict[str, Any] = block.get_mmcif_category("_atom_site.")

        self.arc = ARC(self, atom_site_data)
        # self._coords_array = self.extract_positions(atom_site_data)
        self.arc.atoms.coordinates = self.extract_positions(atom_site_data)

        self.ag: AtomGroupType = self._create_mda()

    def extract_positions(self, atom_site_data: Dict[str, Any]) -> NDArray[np.float32]:
        coords_array: NDArray[np.float32] = np.zeros((self.n_atoms, 3), dtype=np.float32)
        coords_array[:, 0] = atom_site_data.get("Cartn_x")
        coords_array[:, 1] = atom_site_data.get("Cartn_y")
        coords_array[:, 2] = atom_site_data.get("Cartn_z")

        return coords_array

    def _create_mda(self) -> AtomGroupType:
        # Create a structured array to ensure unique values for each combination of resname, resid, and chain_id
        assert self.arc is not None, "arc has not been initialized"
        struct_arr = np.rec.fromarrays(  # type: ignore
            [self.arc.residues.resnames, self.arc.residues.resids, self.arc.chains.ids],
            names=str("resnames, resids, chain_ids"),  # type: ignore
        )

        # Use factorize to get the labels and unique values
        resindices, uniques = pd.factorize(struct_arr)  # type: ignore

        resnames, resids, chain_ids = (uniques["resnames"], uniques["resids"], uniques["chain_ids"])

        # Create a new Universe
        uv: UniverseType = mda.Universe.empty(  # type: ignore
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

        uv.atoms.positions = self.arc.atoms.coordinates

        return uv.atoms

    def to_mda(self) -> AtomGroupType:
        return self.ag

    def to_mol(self) -> MolType:
        assert self.arc is not None, "arc has not been initialized"
        obmol = OBMol()  # type: ignore
        obmol.create_mol(
            self.arc,
            self.structure.connections,  # type: ignore
        )
        assert obmol.mol is not None
        return obmol.mol


class TopologyLoader(BaseLoader):
    def __init__(self, *paths: str):
        file_path: str = paths[0]
        super().__init__(file_path)
        universe = mda.Universe(self.file_path)
        self.ag: AtomGroupType = universe.atoms  # type: ignore
        assert self.ag is not None
        if len(paths) > 1:
            self.ag.universe.load_new(paths[1:], format=None, in_memory=False)  # type: ignore

        self.arc = ARC(self, self.ag)  # positions are set when using mda.Universe

    def to_mda(self) -> AtomGroupType:
        return self.ag

    def to_mol(self) -> MolType:
        assert self.arc is not None, "arc has not been initialized"
        obmol = OBMol()  # type: ignore
        obmol.create_mol(
            self.arc,
            self.structure,
        )
        assert obmol.mol is not None
        return obmol.mol

    @classmethod
    def from_mda(cls, ag: AtomGroupType) -> "TopologyLoader":
        top_loader = cls.__new__(cls)
        top_loader.ag = ag.copy()
        top_loader.ag._u = ag.universe.copy()  # type: ignore
        top_loader.structure = None

        top_loader.arc = ARC(top_loader, top_loader.ag)

        return top_loader
