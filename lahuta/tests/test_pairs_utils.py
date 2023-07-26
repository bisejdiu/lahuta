"""
Tests for the array_utils.py module.

The tests are parameterized to test the functions with different input sizes and ratios of
shared elements between the arrays.
"""
import random
from typing import Callable, List, Tuple

import numpy as np
import pytest
from numpy.typing import NDArray

import lahuta.utils.array_utils as au

TestFuncCallable = Callable[[int, float, float], Tuple[NDArray[np.int32], NDArray[np.int32]]]


# pylint: disable=arguments-out-of-order


def unique_pairs(size: int, start: int = 0):
    total_elements = size * 2

    # Generate unique elements
    unique_elements = np.random.choice(range(start, start + total_elements), total_elements, replace=False)

    # Reshape into desired format
    pairs = unique_elements.reshape(-1, 2)

    return pairs


def _generate_test_data(
    size: int, subset_ratio: float = 0.5, extra_ratio: float = 0.1
) -> Tuple[NDArray[np.int32], NDArray[np.int32]]:
    """
    Generate a pair of arrays where the second array is a subset of the first, plus some additional unique pairs.
    """
    if not 0 <= subset_ratio <= 1 or not 0 <= extra_ratio <= 1:
        raise ValueError("subset_ratio and extra_ratio must be a float between 0 and 1.")

    total_size = size
    if subset_ratio > 0:
        subset_size = int(size * subset_ratio)
        extra_size = int(subset_size * extra_ratio)
        total_size += extra_size

    # Generate unique pairs
    all_data = unique_pairs(total_size)

    # Assign the first 'size' pairs to arr1
    arr1 = all_data[:size]

    if subset_ratio > 0:
        # Assign the first 'subset_size' pairs of arr1 to arr2
        arr2_subset = arr1[:subset_size]  # type: ignore

        # Add additional unique pairs to arr2
        arr2_extra = all_data[size : size + extra_size]  # type: ignore

        arr2 = np.concatenate((arr2_subset, arr2_extra), axis=0)
    else:
        # Assign remaining pairs to arr2
        arr2 = all_data[size:]

    return arr1, arr2


def sort_pairs(arr: NDArray[np.int32]) -> NDArray[np.int32]:
    """
    Sort pairs in an array such that the smaller index is in the first column.

    Args:
        arr (NDArray[np.int32]): The array to sort.

    Returns:
        NDArray[np.int32]: The sorted array.
    """
    # Ensure the smaller index is in the first column
    arr = np.sort(arr, axis=1)

    # Use lexsort to get sorted indices from large to small,
    # then use it to index into the sorted array
    indices = np.lexsort((arr[:, 1], arr[:, 0]))  # type: ignore

    return arr[indices]


def check_mask(mask: NDArray[np.bool_], arr1: NDArray[np.int32], arr2: NDArray[np.int32]) -> None:
    """
    Check the properties of a boolean mask and the associated arrays.
    """
    assert mask.dtype == np.bool_
    assert mask.sum() == np.isin(arr1, arr2).all(axis=1).sum()
    assert np.isin(arr1, arr2).all()


def check_edge_cases(arr1: NDArray[np.int32], arr2: NDArray[np.int32]) -> None:
    """
    Check edge cases of the boolean mask.
    """
    # Test with arrays that have no common elements
    arr1 = arr1 * -1
    mask = au.find_shared_pairs(arr1, arr2)
    assert mask.sum() == 0

    # Test with identical arrays
    arr1 = np.random.randint(0, 10000, size=(1000, 2))  # type: ignore
    arr2 = arr1.copy()
    mask = au.find_shared_pairs(arr1, arr2)
    assert mask.all()


def check_union(
    union_arr: NDArray[np.int32], indices: NDArray[np.int32], arr1: NDArray[np.int32], arr2: NDArray[np.int32]
) -> None:
    """
    Check the properties of a union of arrays and the associated indices.
    """

    unique_concat = np.unique(np.concatenate((arr1, arr2), axis=0), axis=0)

    assert union_arr.shape[0] == unique_concat.shape[0]
    assert indices.shape[0] == unique_concat.shape[0]
    assert np.isin(union_arr, unique_concat).all()
    assert (sort_pairs(union_arr) == sort_pairs(unique_concat)).all()


def check_symmetric_difference(
    mask_a: NDArray[np.bool_], mask_b: NDArray[np.bool_], arr1: NDArray[np.int32], arr2: NDArray[np.int32]
) -> None:
    """
    Check the properties of a symmetric difference of arrays.
    """
    # The masked arrays should not have any common elements
    assert not np.isin(arr1[mask_a], arr2[mask_b]).any()

    # The masked arrays should be equal to the difference in both directions
    assert np.array_equal(arr1[mask_a], arr1[au.difference(arr1, arr2)])
    assert np.array_equal(arr2[mask_b], arr2[au.difference(arr2, arr1)])


@pytest.fixture(scope="module", name="call_func")
def generate_test_data():
    """
    Fixture that returns a function that generates test data.

    Returns:
        function: A function that generates test data.
    """
    return _generate_test_data


# define the parameters for the test
# params = [(1000, 0.5, 0.1), (2000, 0.7, 0.2), (500, 0.3, 0.05)]
params: List[Tuple[int, float, float]] = []
for _ in range(10):
    a = random.randint(500, 20000)  # parameter a: integer between 500 and 2000
    b = round(random.uniform(0.05, 0.95), 2)  # parameter b: float between 0.3 and 0.7, rounded to two decimal places
    c = round(random.uniform(0.05, 0.95), 2)  # parameter c: float between 0.05 and 0.2, rounded to two decimal places
    params.append((a, b, c))


@pytest.mark.parametrize("size,subset_ratio,extra_ratio", params)
def test_shared_pairs(
    call_func: TestFuncCallable,
    size: int,
    subset_ratio: float,
    extra_ratio: float,
) -> None:
    """
    Test the shared pairs function by verifying that the returned boolean mask identifies the correct
    elements in `arr1` that also appear in `arr2`.
    """
    arr1, arr2 = call_func(size, subset_ratio, extra_ratio)
    mask = au.find_shared_pairs(arr1, arr2)
    arr1 = arr1[mask]

    check_mask(mask, arr1, arr2)
    check_edge_cases(arr1, arr2)


@pytest.mark.parametrize("size,subset_ratio,extra_ratio", params)
def test_intersection(
    call_func: TestFuncCallable,
    size: int,
    subset_ratio: float,
    extra_ratio: float,
) -> None:
    """
    Test the intersection function by verifying that the returned boolean mask identifies the correct
    elements in `arr1` that also appear in `arr2`.
    """
    arr1, arr2 = call_func(size, subset_ratio, extra_ratio)
    mask = au.intersection(arr1, arr2)
    arr1 = arr1[mask]

    check_mask(mask, arr1, arr2)
    check_edge_cases(arr1, arr2)


@pytest.mark.parametrize("size,subset_ratio,extra_ratio", params)
def test_difference(
    call_func: TestFuncCallable,
    size: int,
    subset_ratio: float,
    extra_ratio: float,
) -> None:
    """
    Test the difference function by verifying that the returned boolean mask identifies the correct
    elements in `arr1` that do not appear in `arr2`.
    """
    arr1, arr2 = call_func(size, subset_ratio, extra_ratio)

    # Calculate difference
    mask = au.difference(arr1, arr2)

    arr1 = arr1[np.invert(mask)]  # invert mask to get the elements in arr1 that are not in arr2

    check_mask(np.invert(mask), arr1, arr2)
    check_edge_cases(arr1, arr2)


@pytest.mark.parametrize("size,subset_ratio,extra_ratio", params)
def test_union(
    call_func: TestFuncCallable,
    size: int,
    subset_ratio: float,
    extra_ratio: float,
) -> None:
    """
    Test the union function by verifying that the returned arrays contain the correct
    elements from `arr1` and `arr2` and the associated indices.
    """
    arr1, arr2 = call_func(size, subset_ratio, extra_ratio)
    arr1 = sort_pairs(arr1)
    arr2 = sort_pairs(arr2)

    # Calculate union
    union_arr, indices = au.union(arr1, arr2)

    check_union(union_arr, indices, arr1, arr2)


@pytest.mark.parametrize("size,subset_ratio,extra_ratio", params)
def test_symmetric_difference(
    call_func: TestFuncCallable,
    size: int,
    subset_ratio: float,
    extra_ratio: float,
) -> None:
    """
    Test the symmetric difference function by verifying that the returned boolean masks identify the correct
    elements in `arr1` and `arr2`.
    """
    arr1, arr2 = call_func(size, subset_ratio, extra_ratio)
    arr1 = sort_pairs(arr1)
    arr2 = sort_pairs(arr2)

    # Calculate symmetric difference
    mask_a, mask_b = au.symmetric_difference(arr1, arr2)

    check_symmetric_difference(mask_a, mask_b, arr1, arr2)


@pytest.mark.parametrize("size,subset_ratio,extra_ratio", params)
def test_isdisjoint(
    call_func: TestFuncCallable,
    size: int,
    subset_ratio: float,
    extra_ratio: float,
) -> None:
    """
    Test the isdisjoint function by verifying that it returns the correct boolean output
    when given two arrays with and without overlapping elements.
    """
    # Generate disjoint data
    arr1, arr2 = call_func(size, 0, extra_ratio)

    arr1 = sort_pairs(arr1)
    arr2 = sort_pairs(arr2)

    # Check if they are disjoint
    assert au.isdisjoint(arr1, arr2)

    # Generate non-disjoint data
    arr1, arr2 = call_func(size, subset_ratio, extra_ratio)

    arr1 = sort_pairs(arr1)
    arr2 = sort_pairs(arr2)

    # Check if they are not disjoint
    assert not au.isdisjoint(arr1, arr2)


@pytest.mark.parametrize("size,subset_ratio,extra_ratio", params)
def test_issubset(
    call_func: TestFuncCallable,
    size: int,
    subset_ratio: float,
    extra_ratio: float,
) -> None:
    """
    Test the issubset function by verifying that it returns the correct boolean output
    when given an array and a larger array that contains it.
    """
    arr1, arr2 = call_func(size, subset_ratio, 0)
    arr1 = sort_pairs(arr1)
    arr2 = sort_pairs(arr2)

    # Check if arr2 is a subset of arr1
    assert au.issubset(arr2, arr1)

    arr1, arr2 = call_func(size, subset_ratio, extra_ratio)
    arr1 = sort_pairs(arr1)
    arr2 = sort_pairs(arr2)

    # Check if arr2 is not a subset of arr1
    assert not au.issubset(arr2, arr1)


@pytest.mark.parametrize("size,subset_ratio,extra_ratio", params)
def test_issuperset(
    call_func: TestFuncCallable,
    size: int,
    subset_ratio: float,
    extra_ratio: float,
) -> None:
    """
    Test the issuperset function by verifying that it returns the correct boolean output
    when given an array and a smaller array that it contains.
    """
    # Generate test data where arr1 is a superset of arr2
    arr1, arr2 = call_func(size, subset_ratio, 0)

    arr1 = sort_pairs(arr1)
    arr2 = sort_pairs(arr2)

    # Check if arr1 is a superset of arr2
    assert au.issuperset(arr1, arr2)

    # Generate test data where arr1 is not a superset of arr2
    arr1, arr2 = call_func(size, subset_ratio, extra_ratio)

    arr1 = sort_pairs(arr1)
    arr2 = sort_pairs(arr2)

    # Check if arr1 is not a superset of arr2
    assert not au.issuperset(arr1, arr2)


@pytest.mark.parametrize("size,subset_ratio,extra_ratio", params)
def test_isequal(
    call_func: TestFuncCallable,
    size: int,
    subset_ratio: float,
    extra_ratio: float,
) -> None:
    """
    Test the isequal function by verifying that it returns the correct boolean output
    when given two identical arrays and two different arrays.
    """
    # Generate test data where arr1 is equal to arr2
    arr1, arr2 = call_func(size, 1, 0)

    arr1 = sort_pairs(arr1)
    arr2 = sort_pairs(arr2)

    # Check if arr1 is equal to arr2
    assert au.isequal(arr1, arr2)

    # Generate test data where arr1 is not equal to arr2
    arr1, arr2 = call_func(size, subset_ratio, extra_ratio)

    arr1 = sort_pairs(arr1)
    arr2 = sort_pairs(arr2)

    # Check if arr1 is not equal to arr2
    assert not au.isequal(arr1, arr2)


@pytest.mark.parametrize("size,subset_ratio,_", params)
def test_isunique(
    call_func: TestFuncCallable,
    size: int,
    subset_ratio: float,
    _: float,
) -> None:
    """
    Test the isunique function by verifying that it returns the correct boolean output
    when given an array with no duplicates and an array with duplicates.
    """
    # Generate test data with no duplicates
    # arr1 = np.random.randint(0, 10000, size=(1000, 2))
    arr1, arr2 = call_func(size, subset_ratio, 0)

    arr1 = sort_pairs(arr1)
    arr2 = sort_pairs(arr2)

    # Check if arr1 has no duplicates
    assert au.isunique(arr1)
    assert au.isunique(arr2)
    # Check if arr2 has duplicates
    arr3 = np.concatenate((arr1, arr2), axis=0)
    assert not au.isunique(arr3)


@pytest.mark.parametrize("size,subset_ratio,_", params)
def test_is_strict_subset(
    call_func: TestFuncCallable,
    size: int,
    subset_ratio: float,
    _: float,
) -> None:
    """
    Test the is_strict_subset function by verifying that it returns the correct boolean output
    when given an array and a larger array that contains it and when given two identical arrays.
    """

    arr1, arr2 = call_func(size, subset_ratio, 0)
    arr1 = sort_pairs(arr1)
    arr2 = sort_pairs(arr2)

    # Check if arr1 is a strict subset of arr2
    assert au.is_strict_subset(arr2, arr1)

    arr1, arr2 = call_func(size, 1, 0)
    arr1 = sort_pairs(arr1)
    arr2 = sort_pairs(arr2)

    # Check if arr1 is not a strict subset of arr2
    assert not au.is_strict_subset(arr2, arr1)


@pytest.mark.parametrize("size,subset_ratio,_", params)
def test_is_strict_superset(
    call_func: TestFuncCallable,
    size: int,
    subset_ratio: float,
    _: float,
) -> None:
    """
    Test the is_strict_superset function by verifying that it returns the correct boolean output
    when given an array and a smaller array that it contains and when given two identical arrays.
    """
    # Generate test data where arr1 is a strict superset of arr2
    arr1, arr2 = call_func(size, subset_ratio, 0)

    arr1 = sort_pairs(arr1)
    arr2 = sort_pairs(arr2)

    # Check if arr1 is a strict superset of arr2
    assert au.is_strict_superset(arr1, arr2)

    # Generate test data where arr1 is equal to arr2
    arr1, arr2 = call_func(size, 1, 0)

    arr1 = sort_pairs(arr1)
    arr2 = sort_pairs(arr2)

    # Check if arr1 is not a strict superset of arr2
    assert not au.is_strict_superset(arr1, arr2)
