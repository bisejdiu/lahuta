import json
import random
import warnings
from pathlib import Path
from typing import Any, Callable, List, Tuple, Type

import numpy as np
import pytest
from numpy.typing import NDArray

import lahuta.utils.array_utils as au
from lahuta.core.neighbors import NeighborPairs
from lahuta.core.universe import Universe

TestFuncCallable = Callable[[NDArray[np.int32], float, float], NDArray[np.int32]]


def unique_pairs(pairs: NDArray[np.int32], size: int, start: int = 0):
    total_elements = size * 2

    # Generate unique elements excluding existing ones
    upper_bound = start + total_elements + len(np.unique(pairs))
    unique_elements = np.setdiff1d(np.arange(start, upper_bound), np.unique(pairs), assume_unique=True)

    # Select elements randomly
    unique_elements = np.random.choice(unique_elements, total_elements, replace=False)

    # Reshape into desired format
    new_pairs = unique_elements.reshape(-1, 2)

    return new_pairs


# def _generate_test_data(pairs: NDArray[np.int32], subset_size: int, extra_ratio: float = 0):
#     """
#     Modify the 'pairs' array where a subset of pairs are retained, and some additional unique pairs are added.
#     """
#     if not 0 <= extra_ratio <= 1:
#         raise ValueError("extra_ratio must be a float between 0 and 1.")

#     # Calculate the number of extra pairs to add
#     extra_size = int(subset_size * extra_ratio)

#     # Generate unique pairs
#     extra_pairs = unique_pairs(pairs, extra_size)

#     # Retain a subset of original pairs
#     subset_pairs = pairs[:subset_size]

#     # Concatenate the subset of original pairs with the new unique pairs
#     pairs = np.concatenate((subset_pairs, extra_pairs), axis=0)

#     return pairs


def _generate_test_data(pairs: NDArray[np.int32], subset_ratio: float, extra_ratio: float = 0) -> NDArray[np.int32]:
    """
    Modify the 'pairs' array where a subset of pairs are retained, and some additional unique pairs are added.
    """
    if not 0 <= subset_ratio <= 1 or not 0 <= extra_ratio <= 1:
        raise ValueError("subset_ratio and extra_ratio must be a float between 0 and 1.")

    # Calculate the number of pairs to keep and to add
    subset_size = int(pairs.shape[0] * subset_ratio)
    extra_size = int(subset_size * extra_ratio)

    # Generate unique pairs
    extra_pairs = unique_pairs(pairs, extra_size)

    # Retain a subset of original pairs
    subset_pairs = pairs[:subset_size]

    # Concatenate the subset of original pairs with the new unique pairs
    pairs = np.concatenate((subset_pairs, extra_pairs), axis=0)

    return pairs


@pytest.fixture(scope="module", name="call_func")
def generate_test_data() -> TestFuncCallable:
    """
    Fixture that returns a function that generates test data.

    Returns:
        function: A function that generates test data.
    """
    return _generate_test_data


@pytest.fixture(scope="session")
def data_loader() -> NeighborPairs:
    """
    Fixture to load the data for the tests.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        path_obj = Path(__file__).parent / "data" / "1KX2.pdb"
        universe = Universe(str(path_obj))
    ns = universe.compute_neighbors()
    return ns


# params: List[Tuple[int, float, float]] = []
# for _ in range(2):
#     a = round(random.uniform(0.05, 0.95), 2)
#     b = round(random.uniform(0.05, 0.95), 2)
#     params.append((a, b))


@pytest.fixture(scope="session", name="ns")
def neighbors(data_loader: NeighborPairs) -> NeighborPairs:
    """Helper fixture to get neighbor pairs."""
    return data_loader


# define the parameters for the test
params = [(0.1, 0), (0.2, 0)]


@pytest.mark.parametrize("subset_ratio,extra_ratio", params)
def test_intersection(
    call_func: TestFuncCallable,
    subset_ratio: float,
    extra_ratio: float,
    ns: NeighborPairs,
) -> None:
    """
    Test the intersection function by verifying that the returned boolean mask identifies the correct
    elements in `arr1` that also appear in `arr2`.
    """
    new_pairs = call_func(ns.pairs, subset_ratio, 0)
    new_ns = ns.clone(new_pairs, ns.distances[new_pairs.shape[0]])

    resulting_ns = ns & new_ns
    mask = au.intersection(ns.pairs, new_pairs)
    ref_intersect_pairs = ns.pairs[mask]

    assert (resulting_ns.pairs == new_pairs).all()
    assert (resulting_ns.pairs == ref_intersect_pairs).all()
