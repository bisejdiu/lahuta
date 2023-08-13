import warnings
from pathlib import Path
from typing import Callable, Optional, Tuple

import MDAnalysis as mda
import numpy as np
import pytest
from _pytest.fixtures import FixtureRequest

from lahuta import Luni

# from lahuta.contacts import F
from lahuta.contacts import contacts as C
from lahuta.core.neighbors import NeighborPairs
from lahuta.lahuta_types.mdanalysis import UniverseType

HISTIDINE_RESNAMES = ["HIS", "HID", "HIE", "HIP"]
AROMATIC_RESNAMES = ["PHE", "TYR", "TRP"] + HISTIDINE_RESNAMES

# pylint: disable=attribute-defined-outside-init
# pylint: disable=redefined-outer-name

pytestmark = pytest.mark.ag


class ContactType:
    def __init__(
        self,
        name: str,
        func: Callable[[NeighborPairs], NeighborPairs],
        universe: "UniverseWrapper",
    ) -> None:
        self.name = name
        self.func = func
        self.universe = universe
        self.neighbors_ref: Optional[NeighborPairs] = None
        self.neighbors: Optional[NeighborPairs] = None
        self.neighbors_diff: Optional[NeighborPairs] = None

    def compute_neighbors(self, res_dif: int) -> None:
        self.neighbors_ref = self.universe.u_ref.compute_neighbors(res_dif=res_dif)
        self.neighbors = self.universe.u.compute_neighbors(res_dif=res_dif)

    def compute_diff(self) -> None:
        assert self.neighbors is not None
        assert self.neighbors_ref is not None
        self.neighbors_diff = self.func(self.neighbors_ref) - self.func(self.neighbors)

    def pairs(self) -> Tuple[int, ...]:
        assert self.neighbors is not None
        return self.func(self.neighbors).pairs.shape

    def pairs_ref(self) -> Tuple[int, ...]:
        assert self.neighbors_ref is not None
        return self.func(self.neighbors_ref).pairs.shape


@pytest.fixture(scope="session")
def universe_ref(mda_universe: UniverseType) -> Luni:
    with warnings.catch_warnings(record=True) as _:
        return Luni(mda_universe.atoms)


class UniverseWrapper:
    def __init__(self, mda_u: UniverseType, selection: str, u_ref: Luni) -> None:
        self.mda_u = mda_u
        resnames = self.mda_u.select_atoms(f"all and not ({selection})").residues.resnames
        self.unique_resnames = np.unique(resnames)

        self.u_ref = u_ref
        self.u = Luni(self.mda_u.select_atoms(selection).atoms)


@pytest.fixture(scope="session")
def mda_universe() -> UniverseType:
    pdb_path = Path(__file__).parent / "data" / "1KX2.pdb"
    with warnings.catch_warnings(record=True) as _:
        return mda.Universe(str(pdb_path))  # type: ignore


selections_res_difs = [
    ("all", 1),
    ("protein and not resname ARG", 2),
    ("protein and not resname LYS", 3),
    (f"resname {' '.join(AROMATIC_RESNAMES)} or resname HEC", 2),
    ("resid 1 to 50", 4),
    ("resid 1 to 50 or resname HEC", 2),
    ("(protein and not resname ARG LYS) or resname HEC", 2),
]


class TestMDAnalysis:
    @pytest.fixture(params=selections_res_difs, autouse=True)
    def setup_method(self, request: FixtureRequest, mda_universe: UniverseType, universe_ref: Luni) -> None:
        selection, res_dif = request.param
        with warnings.catch_warnings(record=True) as _:
            self.universe = UniverseWrapper(mda_universe, selection, universe_ref)
            self.mda = self.universe.u_ref.to("mda")
            self.selection = selection

        self.contact_types = [
            ContactType("covalent", C.covalent_neighbors, self.universe),
            ContactType("metalic", C.metalic_neighbors, self.universe),
            ContactType("carbonyl", C.carbonyl_neighbors, self.universe),
            ContactType("hbond", C.hbond_neighbors, self.universe),
            ContactType("weak_hbond", C.weak_hbond_neighbors, self.universe),
            ContactType("ionic", C.ionic_neighbors, self.universe),
            ContactType("aromatic", C.aromatic_neighbors, self.universe),
            ContactType("hydrophobic", C.hydrophobic_neighbors, self.universe),
            ContactType("polar_hbond", C.polar_hbond_neighbors, self.universe),
            ContactType("polar_weak_hbond", C.weak_polar_hbond_neighbors, self.universe),
            ContactType("vdw", C.vdw_neighbors, self.universe),
            # ContactType("plane_plane", F.plane_plane_neighbors, self.universe),
        ]

        for contact_type in self.contact_types:
            contact_type.compute_neighbors(res_dif)
            contact_type.compute_diff()

    def test_subset_check(self) -> None:
        for contact_type in self.contact_types:
            message = f"{contact_type.name} neighbors are not a subset of the reference neighbors"
            assert contact_type.neighbors is not None, "Neighbors are None"
            assert contact_type.neighbors_ref is not None
            assert contact_type.neighbors.issubset(contact_type.neighbors_ref), message

    def test_number_of_pairs(self) -> None:
        for contact_type in self.contact_types:
            assert (
                contact_type.pairs() <= contact_type.pairs_ref()
            ), f"{contact_type.name} neighbors pairs are not less than or equal to the reference neighbors pairs"

    def test_correct_residue_content(self) -> None:
        atoms = self.mda.universe.atoms
        assert atoms is not None
        for contact_type in self.contact_types:
            if contact_type.name == "covalent":
                continue
            assert contact_type.neighbors_diff is not None
            for pair in contact_type.neighbors_diff.pairs:
                x1, x2 = atoms[pair[0]], atoms[pair[1]]
                assert (x1.resname in self.universe.unique_resnames) or (
                    x2.resname in self.universe.unique_resnames
                ), f"Pair {pair} with resnames {x1.resname} and {x2.resname} got wrongly picked up by {contact_type.name} neighbors, for selection '{self.selection}', and resnames {self.universe.unique_resnames}"
            assert True
