from pathlib import Path
from typing import List, Set

import MDAnalysis as mda
import numpy as np
import pytest

from lahuta.contacts import contacts as C
from lahuta.core.universe import Universe


class ContactType:
    def __init__(self, name: str, func, universe):
        self.name = name
        self.func = func
        self.universe = universe
        self.neighbors_ref = None
        self.neighbors = None
        self.neighbors_diff = None

    def compute_neighbors(self, res_dif: int):
        self.neighbors_ref = self.universe.u_ref.compute_neighbors(res_dif=res_dif)
        self.neighbors = self.universe.u.compute_neighbors(res_dif=res_dif)

    def compute_diff(self):
        self.neighbors_diff = self.func(self.neighbors_ref) - self.func(self.neighbors)

    def pairs(self):
        return self.func(self.neighbors).pairs.shape

    def pairs_ref(self):
        return self.func(self.neighbors_ref).pairs.shape


class UniverseWrapper:
    def __init__(self, pdb_path: str, selection: str):
        self.mda_u = mda.Universe(pdb_path)
        resnames = self.mda_u.select_atoms(
            f"all and not ({selection})"
        ).residues.resnames
        self.unique_resnames = np.unique(resnames)

        # Initialize *my* Universe class
        self.u_ref = Universe(self.mda_u.atoms)
        self.u = Universe(self.mda_u.select_atoms(selection).atoms)


selections_res_difs = [
    ("protein and not resname ARG", 2),
    # ("protein and not resname LYS", 3),
]


class TestMDAnalysis:
    @pytest.fixture(params=selections_res_difs, autouse=True)
    def setup_method(self, request):
        print("Running setup_method")
        selection, res_dif = request.param
        path_obj = Path(__file__).parent / "data" / "1KX2.pdb"
        self.universe = UniverseWrapper(str(path_obj), selection)
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
            ContactType(
                "polar_weak_hbond", C.weak_polar_hbond_neighbors, self.universe
            ),
            ContactType("vdw", C.vdw_neighbors, self.universe),
        ]

        for contact_type in self.contact_types:
            contact_type.compute_neighbors(res_dif)
            contact_type.compute_diff()

    def test_subset_check(self):
        for contact_type in self.contact_types:
            print("contact_type", contact_type, type(contact_type))
            assert contact_type.neighbors.issubset(
                contact_type.neighbors_ref
            ), f"{contact_type.name} neighbors are not a subset of the reference neighbors"

    def test_number_of_pairs(self):
        for contact_type in self.contact_types:
            assert (
                contact_type.pairs() <= contact_type.pairs_ref()
            ), f"{contact_type.name} neighbors pairs are not less than or equal to the reference neighbors pairs"
