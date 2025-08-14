from itertools import combinations

import numpy as np
import pytest

from lahuta import Luni
from lahuta.core.neighbors import NeighborPairs
from tests.test_data_loading import TEST_PARAMS

DSSP_SELECTION_KW = ["helix_3_10", "pi_helix", "turn", "bend"]  # , "no_ss", "helix"]
DEFAULT_SELECTION_STRINGS = [
    # "all",
    # "protein",
    "backbone",
    "resname ALA VAL ILE LEU",
    "resname ARG LYS ASP GLU and not backbone",
    "hydrophobic and name CA CB",
    "cyclic and backbone",
    "acidic and not backbone",
    "water",
]
SELECTION_STRINGS = DEFAULT_SELECTION_STRINGS + DSSP_SELECTION_KW

SELECTION_STRING_COMBINATIONS = list(combinations(SELECTION_STRINGS, 2))


@pytest.fixture(scope="session", params=TEST_PARAMS)
def data_loader(request: pytest.FixtureRequest) -> tuple[Luni, NeighborPairs]:
    _, pdb_class, pdb_flag = request.param
    luni = Luni(pdb_class(pdb_flag).file_loc)
    return luni, luni.neighbors()


def filter_luni(luni: Luni, target_spec: Luni) -> NeighborPairs:
    return luni.neighbors(target_spec=target_spec) if luni.n_atoms != 0 else NeighborPairs(luni)


@pytest.mark.parametrize("sel_pair", SELECTION_STRING_COMBINATIONS)
def test_luni_filtering(sel_pair: tuple[str, ...], data_loader: tuple[Luni, NeighborPairs]) -> None:
    luni, ns = data_loader

    sel1_luni, sel2_luni = luni.filter(sel_pair[0]), luni.filter(sel_pair[1])
    filtered_ns = filter_luni(sel1_luni, sel2_luni)
    # assert filtered_ns.issubset(ns), f"Subset check failed for {sel_pair}"

    matching_indices = filtered_ns.indices
    assert np.isin(matching_indices, ns.indices).all(), f"Selection pair {sel_pair} failed"
    assert np.isin(matching_indices, sel2_luni.indices).any(axis=1).all(), f"Selection pair {sel_pair} failed"

    reverse_result = filter_luni(sel2_luni, sel1_luni)
    assert np.array_equal(matching_indices, reverse_result.indices), f"Selection pair {sel_pair} failed"
