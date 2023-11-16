import random
from typing import Callable

import pytest

import numpy as np

import lahuta.utils.set_defs as sd
from lahuta.core.neighbors import LabeledNeighborPairs
from tests.test_nb_setdefs import setup_data

T = tuple[LabeledNeighborPairs, LabeledNeighborPairs]

EXPECTED_RESULTS = {
    's1': [
        [('A', '76', 'TYR', 'CZ'), ('A', '154', 'LYS', 'CB')],
        [('A', '40', 'MET', 'SD'), ('A', '45', 'VAL', 'CB')],
        [('A', '139', 'TYR', 'CA'), ('A', '142', 'ILE', 'CG2')],
        [('A', '128', 'ILE', 'CD1'), ('A', '218', 'PRO', 'CA')],
        [('A', '122', 'LEU', 'N'), ('A', '172', 'SER', 'CB')],
        [('A', '86', 'LEU', 'O'), ('A', '89', 'GLY', 'C')],
        [('A', '74', 'THR', 'C'), ('A', '78', 'ILE', 'CD1')],
        [('A', '316', 'PHE', 'CD2'), ('A', '339', 'ASN', 'CA')],
        [('A', '67', 'PHE', 'CZ'), ('A', '363', 'ALA', 'O')],
        [('A', '349', 'ASN', 'N'), ('A', '352', 'ILE', 'CD1')],
    ],
    's2': [
        [('A', '84', 'SER', 'OG'), ('A', '127', 'SER', 'CB')],
        [('A', '47', 'GLY', 'O'), ('A', '50', 'VAL', 'CG1')],
        [('A', '154', 'LYS', 'CA'), ('A', '157', 'SER', 'N')],
        [('A', '136', 'LEU', 'CG'), ('A', '225', 'VAL', 'CG2')],
        [('A', '129', 'TRP', 'CD1'), ('A', '168', 'SER', 'OG')],
        [('A', '95', 'LEU', 'CD1'), ('A', '112', 'LEU', 'CD2')],
        [('A', '81', 'MET', 'CG'), ('A', '349', 'ASN', 'ND2')],
        [('A', '76', 'ASN', 'ND2'), ('A', '154', 'LYS', 'CD')],
        [('A', '296', 'GLN', 'C'), ('A', '300', 'LEU', 'N')],
        [('A', '50', 'VAL', 'CG1'), ('A', '93', 'MET', 'SD')],
    ],
}


@pytest.fixture
def non_transformed_maps(setup_data: T) -> T:
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
            s.partner1
        with pytest.raises(NotImplementedError):
            s.partner2
        with pytest.raises(NotImplementedError):
            s.distances
        with pytest.raises(NotImplementedError):
            s.indices

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

    assert s1.pairs.tolist() == EXPECTED_RESULTS['s1']
    assert s2.pairs.tolist() == EXPECTED_RESULTS['s2']


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
        for subsample_size in range(1, len(unique_values) + 1):
            subsample = np.random.choice(unique_values, subsample_size, replace=False)
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
