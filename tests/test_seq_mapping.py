import random
from typing import Callable

import numpy as np
import pytest

import lahuta.utils.set_defs as sd
from lahuta.core.neighbors import LabeledNeighborPairs
from tests.test_nb_setdefs import setup_data  # noqa: F401

T = tuple[LabeledNeighborPairs, LabeledNeighborPairs]

EXPECTED_RESULTS = {
    "s1": [
        [("A", "CZ", "76", "TYR"), ("A", "CB", "154", "LYS")],
        [("A", "SD", "40", "MET"), ("A", "CB", "45", "VAL")],
        [("A", "CA", "139", "TYR"), ("A", "CG2", "142", "ILE")],
        [("A", "CD1", "128", "ILE"), ("A", "CA", "218", "PRO")],
        [("A", "N", "122", "LEU"), ("A", "CB", "172", "SER")],
        [("A", "O", "86", "LEU"), ("A", "C", "89", "GLY")],
        [("A", "C", "74", "THR"), ("A", "CD1", "78", "ILE")],
        [("A", "CD2", "316", "PHE"), ("A", "CA", "339", "ASN")],
        [("A", "CZ", "67", "PHE"), ("A", "O", "363", "ALA")],
        [("A", "N", "349", "ASN"), ("A", "CD1", "352", "ILE")],
    ],
    "s2": [
        [("A", "OG", "84", "SER"), ("A", "CB", "127", "SER")],
        [("A", "O", "47", "GLY"), ("A", "CG1", "50", "VAL")],
        [("A", "CA", "154", "LYS"), ("A", "N", "157", "SER")],
        [("A", "CG", "136", "LEU"), ("A", "CG2", "225", "VAL")],
        [("A", "CD1", "129", "TRP"), ("A", "OG", "168", "SER")],
        [("A", "CD1", "95", "LEU"), ("A", "CD2", "112", "LEU")],
        [("A", "CG", "81", "MET"), ("A", "ND2", "349", "ASN")],
        [("A", "ND2", "76", "ASN"), ("A", "CD", "154", "LYS")],
        [("A", "C", "296", "GLN"), ("A", "N", "300", "LEU")],
        [("A", "CG1", "50", "VAL"), ("A", "SD", "93", "MET")],
    ],
}


@pytest.fixture()
def non_transformed_maps(setup_data: T) -> T:  # noqa: F811
    s1, s2 = setup_data
    return s1, s2


def run_common_tests(s1: LabeledNeighborPairs, s2: LabeledNeighborPairs) -> None:
    assert sd.check_union(s1, s2)
    assert sd.check_intersection(s1, s2)
    assert sd.check_difference(s1, s2)
    assert sd.check_symmetric_difference(s1, s2)
    assert sd.check_subset(s1)
    assert sd.check_superset(s1)
    assert sd.check_proper_subset(s1)
    assert sd.check_proper_superset(s1)
    assert sd.check_disjoint(s1, s2)


def get_random_subsample_from_s1(s1: LabeledNeighborPairs, size: int = 10) -> LabeledNeighborPairs:
    random.seed(42)
    indices = random.sample(range(s1.pairs.shape[0]), size)
    return LabeledNeighborPairs(s1.pairs[indices])


def get_random_subsample_from_s2(s2: LabeledNeighborPairs, size: int = 10) -> LabeledNeighborPairs:
    random.seed(42)
    indices = random.sample(range(s2.pairs.shape[0]), size)
    return LabeledNeighborPairs(s2.pairs[indices])


# Test functions
def test_size(non_transformed_maps: T) -> None:
    s1, s2 = non_transformed_maps

    assert len(s1) == 10160
    assert s1.pairs.shape == (10160, 2)
    assert len(s2) == 8863
    assert s2.pairs.shape == (8863, 2)


def test_error_checks(non_transformed_maps: T) -> None:
    s1, s2 = non_transformed_maps

    def error_checks(s: LabeledNeighborPairs) -> None:
        with pytest.raises(NotImplementedError):
            s.type_filter(1, 2, 3, 4)
        with pytest.raises(NotImplementedError):
            s.index_filter(1, 2, 3, 4)
        with pytest.raises(NotImplementedError):
            s.distance_filter(1, 2, 3, 4)
        with pytest.raises(NotImplementedError):
            s.numeric_filter(1, 2, 3, 4)
        with pytest.raises(NotImplementedError):
            s.radius_filter(1, 2, 3, 4)
        with pytest.raises(NotImplementedError):
            s.hbond_distance_filter(1, 2, 3, 4)
        with pytest.raises(NotImplementedError):
            s.hbond_angle_filter(1, 2, 3, 4)
        with pytest.raises(NotImplementedError):
            _ = s.partner1
        with pytest.raises(NotImplementedError):
            _ = s.partner2
        with pytest.raises(NotImplementedError):
            _ = s.distances
        with pytest.raises(NotImplementedError):
            _ = s.indices

    error_checks(s1)
    error_checks(s2)


def test_contents(non_transformed_maps: T) -> None:
    s1, s2 = non_transformed_maps

    s1 = get_random_subsample_from_s1(s1)
    s2 = get_random_subsample_from_s2(s2)

    assert s1.pairs.shape == (10, 2)
    assert s2.pairs.shape == (10, 2)

    assert s1 == get_random_subsample_from_s1(s1)
    assert s2 == get_random_subsample_from_s2(s2)

    assert s1 != s2
    assert s2 != s1

    assert s1 != s2.inverse()
    assert s2 != s1.inverse()

    assert s1.pairs.tolist() == EXPECTED_RESULTS["s1"]
    assert s2.pairs.tolist() == EXPECTED_RESULTS["s2"]


class TestHelper:
    __test__ = False  # Utility class for testing. Not a test itself.

    def __init__(self, s: LabeledNeighborPairs, field_names: tuple[str, ...]) -> None:
        self.s = s
        self.field_names = field_names

    def process_elements(self, select_strategy: Callable[..., None]) -> None:
        for select_field in self.field_names:  # chain_ids, resids, resnames, names
            for remove_field in self.field_names:
                if select_field != remove_field:
                    # select all unique values from the select_field
                    unique_values = np.unique(self.s.pairs[select_field])
                    select_strategy(self, select_field, remove_field, unique_values)

    def standard_select(self, select_field: str, remove_field: str, unique_values: list[str]) -> None:
        # test one element at a time
        for value in unique_values:
            self._select_and_assert(select_field, remove_field, [value])

    def random_select(self, select_field: str, remove_field: str, unique_values: list[str]) -> None:
        # test random subsamples
        rng = np.random.default_rng()
        for subsample_size in range(1, len(unique_values) + 1):
            subsample = rng.choice(unique_values, subsample_size, replace=False)
            self._select_and_assert(select_field, remove_field, subsample.tolist())

    def _select_and_assert(self, select_field: str, remove_field: str, select_values: list[str]) -> None:
        selected_s = self.s.select(**{select_field: select_values})
        selected_s = selected_s.remove(remove_field)
        create_mask = np.isin(selected_s.pairs[select_field], select_values)
        assert np.all(
            selected_s.pairs[create_mask][remove_field] == ""
        ), f"Test failed for combination: {select_field}, {remove_field}"


@pytest.mark.parametrize("sample_index", [0, 1])
@pytest.mark.parametrize("sample_size", [1, 5, 10, 25, 50])
def test_single_element_removal(sample_index: int, sample_size: int, non_transformed_maps: T) -> None:
    s = get_random_subsample_from_s1(non_transformed_maps[sample_index], sample_size)
    field_names = s.pairs.dtype.names
    assert field_names is not None

    helper = TestHelper(s, field_names)
    helper.process_elements(TestHelper.standard_select)


@pytest.mark.parametrize("sample_index", [0, 1])
@pytest.mark.parametrize("sample_size", [1, 5, 10, 25, 50])
def test_variable_element_removal(sample_index: int, sample_size: int, non_transformed_maps: T) -> None:
    s = get_random_subsample_from_s1(non_transformed_maps[sample_index], sample_size)
    field_names = s.pairs.dtype.names
    assert field_names is not None

    helper = TestHelper(s, field_names)
    helper.process_elements(TestHelper.random_select)
