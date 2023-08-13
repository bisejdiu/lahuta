import random
import warnings
from pathlib import Path
from typing import Callable, List, Tuple

import numpy as np
import pytest
from numpy.typing import NDArray

import lahuta.utils.array_utils as au
from lahuta.core.luni import Luni
from lahuta.core.neighbors import NeighborPairs

pytestmark = pytest.mark.nb

TestFuncCallable = Callable[[NDArray[np.int32], float, float], NDArray[np.int32]]


def unique_pairs(pairs: NDArray[np.int32], size: int, start: int = 0) -> NDArray[np.int32]:
    """
    Generate unique pairs of integers.

    This function generates unique pairs of integers that are not already present in the 'pairs' array.
    Generated pairs are unique in the sense that they do not contain any repeated elements.

    Args:
        pairs (NDArray[np.int32]): An array of pairs of integers.
        size (int): The number of unique pairs to generate.
        start (int): The starting integer from which to generate unique pairs.

    Returns:
        NDArray[np.int32]: An array of unique pairs of integers.
    """

    total_elements = size * 2

    # Generate unique elements excluding existing ones
    upper_bound = start + total_elements + len(np.unique(pairs))
    unique_elements = np.setdiff1d(np.arange(start, upper_bound), np.unique(pairs), assume_unique=True)

    # Select elements randomly
    unique_elements = np.random.choice(unique_elements, total_elements, replace=False)

    # Reshape into desired format
    new_pairs = unique_elements.reshape(-1, 2)

    return new_pairs


def _generate_test_data(pairs: NDArray[np.int32], subset_ratio: float, extra_ratio: float = 0) -> NDArray[np.int32]:
    """
    Generate test data.

    This function generates test data by retaining a subset of the original pairs and adding new unique pairs.
    The way the new pairs are generated is by first generating unique pairs that are not already present in the
    original pairs and then concatenating them with the subset of original pairs.

    Args:
        pairs (NDArray[np.int32]): An array of pairs of integers.
        subset_ratio (float): The ratio of pairs to retain from the original pairs.
        extra_ratio (float): The ratio of extra pairs to add to the original pairs.

    Returns:
        NDArray[np.int32]: An array of pairs of integers.
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
        universe = Luni(str(path_obj))
    ns = universe.compute_neighbors()
    return ns


params: List[Tuple[float, float]] = []
for _ in range(10):
    subset = round(random.uniform(0, 0.5), 2)
    extra = round(random.uniform(0, 0.5), 2)
    params.append((subset, extra))


@pytest.fixture(scope="session", name="ns")
def neighbors(data_loader: NeighborPairs) -> NeighborPairs:
    """Helper fixture to get neighbor pairs."""
    return data_loader


@pytest.mark.parametrize("subset_ratio,_", params)
def test_intersection(
    call_func: TestFuncCallable,
    subset_ratio: float,
    _: float,
    ns: NeighborPairs,
) -> None:
    """
    Test the intersection of two NeighborPairs objects.

    Args:
        call_func (TestFuncCallable): A function that generates test data.
        subset_ratio (float): The ratio of pairs to retain from the original pairs.
        _ (float): Unused.
        ns (NeighborPairs): The NeighborPairs object to test.

    Returns:
        None

    """
    new_pairs = call_func(ns.pairs, subset_ratio, 0)
    new_ns = ns.clone(new_pairs, ns.distances[np.arange(new_pairs.shape[0])])

    resulting_ns = ns & new_ns
    mask = au.intersection(ns.pairs, new_pairs)
    ref_intersect_pairs = ns.pairs[mask]

    assert (resulting_ns.pairs == new_pairs).all()
    assert (resulting_ns.pairs == ref_intersect_pairs).all()


@pytest.mark.parametrize("subset_ratio,extra_ratio", params)
def test_union(
    call_func: TestFuncCallable,
    subset_ratio: float,
    extra_ratio: float,
    ns: NeighborPairs,
) -> None:
    """
    Test the union of two NeighborPairs objects.

    Args:
        call_func (TestFuncCallable): A function that generates test data.
        subset_ratio (float): The ratio of pairs to retain from the original pairs.
        extra_ratio (float): The ratio of extra pairs to add to the original pairs.
        ns (NeighborPairs): The NeighborPairs object to test.

    Returns:
        None
    """
    new_pairs = call_func(ns.pairs, subset_ratio, extra_ratio)
    new_ns = ns.clone(new_pairs, ns.distances[np.arange(new_pairs.shape[0])])

    resulting_ns = ns + new_ns
    ref_union_pairs, _ = au.union(ns.pairs, new_pairs)

    assert (resulting_ns.pairs == ref_union_pairs).all()


@pytest.mark.parametrize("subset_ratio,extra_ratio", params)
def test_difference(
    call_func: TestFuncCallable,
    subset_ratio: float,
    extra_ratio: float,
    ns: NeighborPairs,
) -> None:
    """
    Test the difference of two NeighborPairs objects.

    Args:
        call_func (TestFuncCallable): A function that generates test data.
        size (int): The size of the original pairs.
        subset_ratio (float): The ratio of pairs to retain from the original pairs.
        extra_ratio (float): The ratio of extra pairs to add to the original pairs.
        ns (NeighborPairs): The NeighborPairs object to test.

    Returns:
        None
    """
    new_pairs = call_func(ns.pairs, subset_ratio, extra_ratio)
    new_ns = ns.clone(new_pairs, ns.distances[np.arange(new_pairs.shape[0])])

    resulting_ns = ns - new_ns
    mask = au.difference(ns.pairs, new_pairs)
    ref_diff_pairs = ns.pairs[mask]

    assert (resulting_ns.pairs == ref_diff_pairs).all()


@pytest.mark.parametrize("subset_ratio,extra_ratio", params)
def test_symmetric_difference(
    call_func: TestFuncCallable,
    subset_ratio: float,
    extra_ratio: float,
    ns: NeighborPairs,
) -> None:
    """
    Test the symmetric difference of two NeighborPairs objects.

    Args:
        call_func (TestFuncCallable): A function that generates test data.
        size (int): The size of the original pairs.
        subset_ratio (float): The ratio of pairs to retain from the original pairs.
        extra_ratio (float): The ratio of extra pairs to add to the original pairs.
        ns (NeighborPairs): The NeighborPairs object to test.

    Returns:
        None
    """
    new_pairs = call_func(ns.pairs, subset_ratio, extra_ratio)
    new_ns = ns.clone(new_pairs, ns.distances[np.arange(new_pairs.shape[0])])

    resulting_ns = ns | new_ns
    mask_a, mask_b = au.symmetric_difference(ns.pairs, new_pairs)

    ref_sym_diff_pairs = np.concatenate((ns.pairs[mask_a], new_pairs[mask_b]), axis=0)

    assert (resulting_ns.pairs == ref_sym_diff_pairs).all()


@pytest.mark.parametrize("subset_ratio,extra_ratio", params)
def test_isdisjoint(
    call_func: TestFuncCallable,
    subset_ratio: float,
    extra_ratio: float,
    ns: NeighborPairs,
) -> None:
    """
    Test the isdisjoint method of two NeighborPairs objects.

    Args:
        call_func (TestFuncCallable): A function that generates test data.
        size (int): The size of the original pairs.
        subset_ratio (float): The ratio of pairs to retain from the original pairs.
        extra_ratio (float): The ratio of extra pairs to add to the original pairs.
        ns (NeighborPairs): The NeighborPairs object to test.

    Returns:
        None
    """
    new_pairs = call_func(ns.pairs, subset_ratio, extra_ratio)
    new_ns = ns.clone(new_pairs, ns.distances[np.arange(new_pairs.shape[0])])

    assert ns.isdisjoint(new_ns) == au.isdisjoint(ns.pairs, new_pairs)


@pytest.mark.parametrize("subset_ratio,extra_ratio", params)
def test_issubset(
    call_func: TestFuncCallable,
    subset_ratio: float,
    extra_ratio: float,
    ns: NeighborPairs,
) -> None:
    """
    Test the issubset method of two NeighborPairs objects.

    Args:
        call_func (TestFuncCallable): A function that generates test data.
        size (int): The size of the original pairs.
        subset_ratio (float): The ratio of pairs to retain from the original pairs.
        extra_ratio (float): The ratio of extra pairs to add to the original pairs.
        ns (NeighborPairs): The NeighborPairs object to test.

    Returns:
        None
    """
    new_pairs = call_func(ns.pairs, subset_ratio, extra_ratio)
    new_ns = ns.clone(new_pairs, ns.distances[np.arange(new_pairs.shape[0])])

    assert ns.issubset(new_ns) == au.issubset(ns.pairs, new_pairs)


@pytest.mark.parametrize("subset_ratio,extra_ratio", params)
def test_issuperset(
    call_func: TestFuncCallable,
    subset_ratio: float,
    extra_ratio: float,
    ns: NeighborPairs,
) -> None:
    """
    Test the issuperset method of two NeighborPairs objects.

    Args:
        call_func (TestFuncCallable): A function that generates test data.
        size (int): The size of the original pairs.
        subset_ratio (float): The ratio of pairs to retain from the original pairs.
        extra_ratio (float): The ratio of extra pairs to add to the original pairs.
        ns (NeighborPairs): The NeighborPairs object to test.

    Returns:
        None
    """
    new_pairs = call_func(ns.pairs, subset_ratio, extra_ratio)
    new_ns = ns.clone(new_pairs, ns.distances[np.arange(new_pairs.shape[0])])

    assert ns.issuperset(new_ns) == au.issuperset(ns.pairs, new_pairs)


@pytest.mark.parametrize("subset_ratio,extra_ratio", params)
def test_isequal(
    call_func: TestFuncCallable,
    subset_ratio: float,
    extra_ratio: float,
    ns: NeighborPairs,
) -> None:
    """
    Test the isequal method of two NeighborPairs objects.

    Args:
        call_func (TestFuncCallable): A function that generates test data.
        size (int): The size of the original pairs.
        subset_ratio (float): The ratio of pairs to retain from the original pairs.
        extra_ratio (float): The ratio of extra pairs to add to the original pairs.
        ns (NeighborPairs): The NeighborPairs object to test.

    Returns:
        None
    """
    new_pairs = call_func(ns.pairs, subset_ratio, extra_ratio)
    new_ns = ns.clone(new_pairs, ns.distances[np.arange(new_pairs.shape[0])])

    assert ns.isequal(new_ns) == au.isequal(ns.pairs, new_pairs)


@pytest.mark.parametrize("subset_ratio,extra_ratio", params)
def test_isunique(
    call_func: TestFuncCallable,
    subset_ratio: float,
    extra_ratio: float,
    ns: NeighborPairs,
) -> None:
    """
    Test the isunique method of two NeighborPairs objects.

    Args:
        call_func (TestFuncCallable): A function that generates test data.
        size (int): The size of the original pairs.
        subset_ratio (float): The ratio of pairs to retain from the original pairs.
        extra_ratio (float): The ratio of extra pairs to add to the original pairs.
        ns (NeighborPairs): The NeighborPairs object to test.

    Returns:
        None
    """
    new_pairs = call_func(ns.pairs, subset_ratio, extra_ratio)
    new_ns = ns.clone(new_pairs, ns.distances[np.arange(new_pairs.shape[0])])

    # ns.isunique() does not take any arguments
    assert ns.isunique() == au.isunique(ns.pairs)
    assert new_ns.isunique() == au.isunique(new_ns.pairs)


@pytest.mark.parametrize("subset_ratio,extra_ratio", params)
def test_is_strict_subset(
    call_func: TestFuncCallable,
    subset_ratio: float,
    extra_ratio: float,
    ns: NeighborPairs,
) -> None:
    """
    Test the is_strict_subset method of two NeighborPairs objects.

    Args:
        call_func (TestFuncCallable): A function that generates test data.
        subset_ratio (float): The ratio of pairs to retain from the original pairs.
        extra_ratio (float): The ratio of extra pairs to add to the original pairs.
        ns (NeighborPairs): The NeighborPairs object to test.
    """
    new_pairs = call_func(ns.pairs, subset_ratio, extra_ratio)
    new_ns = ns.clone(new_pairs, ns.distances[np.arange(new_pairs.shape[0])])
    assert ns.is_strict_subset(new_ns) == au.is_strict_subset(ns.pairs, new_pairs)


@pytest.mark.parametrize("subset_ratio,extra_ratio", params)
def test_is_strict_superset(
    call_func: TestFuncCallable,
    subset_ratio: float,
    extra_ratio: float,
    ns: NeighborPairs,
) -> None:
    """
    Test the is_strict_superset method of two NeighborPairs objects.

    Args:
        call_func (TestFuncCallable): A function that generates test data.
        subset_ratio (float): The ratio of pairs to retain from the original pairs.
        extra_ratio (float): The ratio of extra pairs to add to the original pairs.
        ns (NeighborPairs): The NeighborPairs object to test.
    """
    new_pairs = call_func(ns.pairs, subset_ratio, extra_ratio)
    new_ns = ns.clone(new_pairs, ns.distances[np.arange(new_pairs.shape[0])])
    assert ns.is_strict_superset(new_ns) == au.is_strict_superset(ns.pairs, new_pairs)


@pytest.mark.parametrize("subset_ratio,_", params)
def test_contains(
    call_func: TestFuncCallable,
    subset_ratio: float,
    _: float,
    ns: NeighborPairs,
) -> None:
    """
    Test the contains method of two NeighborPairs objects.

    Args:
        call_func (TestFuncCallable): A function that generates test data.
        subset_ratio (float): The ratio of pairs to retain from the original pairs.
        _ (float): Unused.
        ns (NeighborPairs): The NeighborPairs object to test.
    """
    new_pairs = call_func(ns.pairs, subset_ratio, 0)
    new_ns = ns.clone(new_pairs, ns.distances[np.arange(new_pairs.shape[0])])
    assert new_ns in ns


@pytest.mark.parametrize("_", params)
def test_eq(
    call_func: TestFuncCallable,
    _: float,
    ns: NeighborPairs,
) -> None:
    """
    Test the eq method of two NeighborPairs objects.

    Args:
        call_func (TestFuncCallable): A function that generates test data.
        _ (float): Unused.
        ns (NeighborPairs): The NeighborPairs object to test.
    """

    new_pairs = call_func(ns.pairs, 1, 0)
    new_ns = ns.clone(new_pairs, ns.distances[np.arange(new_pairs.shape[0])])
    assert new_ns == ns
