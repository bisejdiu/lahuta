import warnings
from pathlib import Path
from typing import Callable, Optional

import MDAnalysis as mda
import numpy as np
import pytest
from _pytest.fixtures import FixtureRequest

from lahuta import Luni

from lahuta.contacts import contacts as C
from lahuta.core.neighbors import NeighborPairs
from lahuta._types.mdanalysis import UniverseType


pytestmark = pytest.mark.trajs

HISTIDINE_RESNAMES = ["HIS", "HID", "HIE", "HIP"]
AROMATIC_RESNAMES = ["PHE", "TYR", "TRP"] + HISTIDINE_RESNAMES


class ContactType:
    """Helper class to compute and store the neighbors of a given type"""

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
        """Compute the neighbors of the given type"""
        self.neighbors_ref = self.universe.u_ref.compute_neighbors(res_dif=res_dif)
        self.neighbors = self.universe.u.compute_neighbors(res_dif=res_dif)

    def compute_diff(self) -> None:
        """Compute the difference between the neighbors of the given type"""
        assert self.neighbors is not None
        assert self.neighbors_ref is not None
        self.neighbors_diff = self.func(self.neighbors_ref) - self.func(self.neighbors)

    def pairs(self) -> tuple[int, ...]:
        """Return the number of pairs"""
        assert self.neighbors is not None
        return self.func(self.neighbors).pairs.shape

    def pairs_ref(self) -> tuple[int, ...]:
        """Return the number of pairs in the reference"""
        assert self.neighbors_ref is not None
        return self.func(self.neighbors_ref).pairs.shape


@pytest.fixture(scope="session")
def universe_ref(mda_universe: UniverseType) -> Luni:
    """Create a reference universe"""
    with warnings.catch_warnings(record=True) as _:
        return Luni(mda_universe.atoms)


class UniverseWrapper:
    """Helper class to store the universe and the selection"""

    def __init__(self, mda_u: UniverseType, selection: str, u_ref: Luni) -> None:
        self.mda_u = mda_u
        resnames = self.mda_u.select_atoms(f"all and not ({selection})").residues.resnames
        self.unique_resnames = np.unique(resnames)

        self.u_ref = u_ref
        self.u = Luni(self.mda_u.select_atoms(selection).atoms)
        self.u.assing_atom_types()


@pytest.fixture(scope="session")
def mda_universe(request: FixtureRequest) -> UniverseType:
    """Create a MDAnalysis universe"""
    use_large_files = request.config.getoption("--large-files")
    if use_large_files:
        coords = Path(__file__).parent / "data" / "197.pdb"
        traj = Path(__file__).parent / "data" / "197.xtc"
    else:
        coords = Path(__file__).parent / "data" / "conf0_197_sel.pdb"
        traj = Path(__file__).parent / "data" / "trj_197_sel.xtc"

    with warnings.catch_warnings(record=True) as _:
        return mda.Universe(str(coords), str(traj))  # type: ignore


selections_res_difs = [
    ("protein or resname TIP3", 2),
    ("protein or resname POPC", 2),
    ("protein or resname POPG", 2),
    ("protein or resname SOD", 2),
    ("(protein and not resname ARG LYS) or resname POPC POPG", 2),
]


class TestMDAnalysis:
    """Test the MDAnalysis implementation of the contacts"""

    @pytest.fixture(params=selections_res_difs, autouse=True)
    def setup_method(self, request: FixtureRequest, mda_universe: UniverseType, universe_ref: Luni) -> None:
        """Setup the test"""
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
        ]

        for contact_type in self.contact_types:
            contact_type.compute_neighbors(res_dif)
            contact_type.compute_diff()

    def test_subset_check(self) -> None:
        """Test that the neighbors are a subset of the reference neighbors"""
        for contact_type in self.contact_types:
            message = f"{contact_type.name} neighbors are not a subset of the reference neighbors"
            assert contact_type.neighbors is not None, "Neighbors are None"
            assert contact_type.neighbors_ref is not None
            assert contact_type.neighbors.issubset(contact_type.neighbors_ref), message

    def test_number_of_pairs(self) -> None:
        """Test that the number of pairs is correct"""
        for contact_type in self.contact_types:
            assert (
                contact_type.pairs() <= contact_type.pairs_ref()
            ), f"{contact_type.name} neighbors pairs are not less than or equal to the reference neighbors pairs"

    def test_correct_residue_content(self) -> None:
        """Test that the pairs are correct"""
        atoms = self.mda.universe.atoms
        assert atoms is not None
        for contact_type in self.contact_types:
            if contact_type.name == "covalent":
                continue
            assert contact_type.neighbors_diff is not None
            for pair in contact_type.neighbors_diff.pairs:
                x1, x2 = atoms[pair[0]], atoms[pair[1]]
                assert (x1.resname in self.universe.unique_resnames) or (x2.resname in self.universe.unique_resnames), (
                    f"Pair {pair} with resnames {x1.resname} and {x2.resname} "
                    f"got wrongly picked up by {contact_type.name} neighbors, "
                    f"for selection '{self.selection}', and resnames {self.universe.unique_resnames}"
                )
            assert True
