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

# pylint: disable=arguments-out-of-order

pytestmark = pytest.mark.au

TestFuncCallable = Callable[[int, float, float], Tuple[NDArray[np.int32], NDArray[np.int32]]]

def unique_pairs(size: int, start: int = 0):
    """
    Generate an array of unique pairs.

    We generate an array of unique pairs by first generating an array of unique elements and then reshaping
    it into an array of pairs. This way we can ensure that the pairs are unique.

    Args:
        size (int): The size of the array to generate.
        start (int): The starting value of the elements.

    Returns:
        NDArray[np.int32]: An array of unique pairs.

    """
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
    Generate test data for the array_utils tests.

    This function generates two arrays of size `size` where `subset_ratio` of the elements in the first array
    are also in the second array. Additionally, `extra_ratio` of the elements in the second array are not in
    the first array.

    Args:
        size (int): The size of the arrays to generate.
        subset_ratio (float): The ratio of elements in `arr1` that are also in `arr2`.
        extra_ratio (float): The ratio of additional elements in `arr2` that are not in `arr1`.

    Returns:
        Tuple[NDArray[np.int32], NDArray[np.int32]]: The two arrays.
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


def check_mask(mask: NDArray[np.bool_], arr1: NDArray[np.int32], arr2: NDArray[np.int32]) -> None:
    """
    Check the properties of the boolean mask.

    Args:
        mask (NDArray[np.bool_]): The boolean mask.
        arr1 (NDArray[np.int32]): The first array.
        arr2 (NDArray[np.int32]): The second array.

    Returns:
        None
    """
    assert mask.dtype == np.bool_
    assert mask.sum() == np.isin(arr1, arr2).all(axis=1).sum()
    assert np.isin(arr1, arr2).all()


def check_edge_cases(arr1: NDArray[np.int32], arr2: NDArray[np.int32]) -> None:
    """
    Check the properties of the arrays when they contain edge cases.

    Args:
        arr1 (NDArray[np.int32]): The first array.
        arr2 (NDArray[np.int32]): The second array.

    Returns:
        None
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
    Check the properties of a union of arrays.

    Args:
        union_arr (NDArray[np.int32]): The union of arr1 and arr2.
        indices (NDArray[np.int32]): The indices of the elements in the union array.
        arr1 (NDArray[np.int32]): The first array.
        arr2 (NDArray[np.int32]): The second array.

    Returns:
        None
    """

    unique_concat = np.unique(np.concatenate((arr1, arr2), axis=0), axis=0)

    assert union_arr.shape[0] == unique_concat.shape[0]
    assert indices.shape[0] == unique_concat.shape[0]
    assert np.isin(union_arr, unique_concat).all()
    assert (au.sort_pairs(union_arr) == au.sort_pairs(unique_concat)).all()


def check_symmetric_difference(
    mask_a: NDArray[np.bool_], mask_b: NDArray[np.bool_], arr1: NDArray[np.int32], arr2: NDArray[np.int32]
) -> None:
    """
    Check the properties of a symmetric difference of arrays.

    Args:
        mask_a (NDArray[np.bool_]): The boolean mask for arr1.
        mask_b (NDArray[np.bool_]): The boolean mask for arr2.
        arr1 (NDArray[np.int32]): The first array.
        arr2 (NDArray[np.int32]): The second array.

    Returns:
        None
    """
    # The masked arrays should not have any common elements
    assert not np.isin(arr1[mask_a], arr2[mask_b]).any()

    # The masked arrays should be equal to the difference in both directions
    assert np.array_equal(arr1[mask_a], arr1[au.difference(arr1, arr2)])
    assert np.array_equal(arr2[mask_b], arr2[au.difference(arr2, arr1)])


@pytest.fixture(scope="module", name="call_func")
def generate_test_data() -> TestFuncCallable:
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
    a = random.randint(500, 20000)
    b = round(random.uniform(0.05, 0.95), 2)
    c = round(random.uniform(0.05, 0.95), 2)
    params.append((a, b, c))


@pytest.mark.parametrize("size,subset_ratio,extra_ratio", params)
def test_shared_pairs(
    call_func: TestFuncCallable,
    size: int,
    subset_ratio: float,
    extra_ratio: float,
) -> None:
    """
    Test the find_shared_pairs function by verifying that the returned boolean mask identifies the correct
    elements in `arr1` that also appear in `arr2`.

    Args:
        call_func (TestFuncCallable): The function that generates test data.
        size (int): The size of the arrays to generate.
        subset_ratio (float): The ratio of elements in `arr1` that are also in `arr2`.
        extra_ratio (float): The ratio of additional elements in `arr2` that are not in `arr1`.

    Returns:
        None
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

    Args:
        call_func (TestFuncCallable): The function that generates test data.
        size (int): The size of the arrays to generate.
        subset_ratio (float): The ratio of elements in `arr1` that are also in `arr2`.
        extra_ratio (float): The ratio of additional elements in `arr2` that are not in `arr1`.

    Returns:
        None
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
    elements in `arr1` that are not in `arr2`.

    Args:
        call_func (TestFuncCallable): The function that generates test data.
        size (int): The size of the arrays to generate.
        subset_ratio (float): The ratio of elements in `arr1` that are also in `arr2`.
        extra_ratio (float): The ratio of additional elements in `arr2` that are not in `arr1`.

    Returns:
        None
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
    Test the union function by verifying that the returned boolean mask identifies the correct
    elements in `arr1` and `arr2`.

    Args:
        call_func (TestFuncCallable): The function that generates test data.
        size (int): The size of the arrays to generate.
        subset_ratio (float): The ratio of elements in `arr1` that are also in `arr2`.
        extra_ratio (float): The ratio of additional elements in `arr2` that are not in `arr1`.

    Returns:
        None
    """
    arr1, arr2 = call_func(size, subset_ratio, extra_ratio)
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
    Test the symmetric_difference function by verifying that the returned boolean masks identify the correct
    elements in `arr1` and `arr2`.

    Args:
        call_func (TestFuncCallable): The function that generates test data.
        size (int): The size of the arrays to generate.
        subset_ratio (float): The ratio of elements in `arr1` that are also in `arr2`.
        extra_ratio (float): The ratio of additional elements in `arr2` that are not in `arr1`.

    Returns:
        None
    """
    arr1, arr2 = call_func(size, subset_ratio, extra_ratio)
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
    when given two disjoint arrays and two non-disjoint arrays.

    Args:
        call_func (TestFuncCallable): The function that generates test data.
        size (int): The size of the arrays to generate.
        subset_ratio (float): The ratio of elements in `arr1` that are also in `arr2`.
        extra_ratio (float): The ratio of additional elements in `arr2` that are not in `arr1`.

    Returns:
        None
    """
    # Generate disjoint data
    arr1, arr2 = call_func(size, 0, extra_ratio)
    # Check if they are disjoint
    assert au.isdisjoint(arr1, arr2)

    # Generate non-disjoint data
    arr1, arr2 = call_func(size, subset_ratio, extra_ratio)
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

    Args:
        call_func (TestFuncCallable): The function that generates test data.
        size (int): The size of the arrays to generate.
        subset_ratio (float): The ratio of elements in `arr1` that are also in `arr2`.
        extra_ratio (float): The ratio of additional elements in `arr2` that are not in `arr1`.

    Returns:
        None
    """
    arr1, arr2 = call_func(size, subset_ratio, 0)
    # Check if arr2 is a subset of arr1
    assert au.issubset(arr2, arr1)

    arr1, arr2 = call_func(size, subset_ratio, extra_ratio)
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

    Args:
        call_func (TestFuncCallable): The function that generates test data.
        size (int): The size of the arrays to generate.
        subset_ratio (float): The ratio of elements in `arr1` that are also in `arr2`.
        extra_ratio (float): The ratio of additional elements in `arr2` that are not in `arr1`.

    Returns:
        None
    """
    # Generate test data where arr1 is a superset of arr2
    arr1, arr2 = call_func(size, subset_ratio, 0)
    # Check if arr1 is a superset of arr2
    assert au.issuperset(arr1, arr2)

    # Generate test data where arr1 is not a superset of arr2
    arr1, arr2 = call_func(size, subset_ratio, extra_ratio)
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
    when given two identical arrays and two non-identical arrays.

    Args:
        call_func (TestFuncCallable): The function that generates test data.
        size (int): The size of the arrays to generate.
        subset_ratio (float): The ratio of elements in `arr1` that are also in `arr2`.
        extra_ratio (float): The ratio of additional elements in `arr2` that are not in `arr1`.

    Returns:
        None
    """
    # Generate test data where arr1 is equal to arr2
    arr1, arr2 = call_func(size, 1, 0)
    assert au.isequal(arr1, arr2)

    # Generate test data where arr1 is not equal to arr2
    arr1, arr2 = call_func(size, subset_ratio, extra_ratio)
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

    Args:
        call_func (TestFuncCallable): The function that generates test data.
        size (int): The size of the arrays to generate.
        subset_ratio (float): The ratio of elements in `arr1` that are also in `arr2`.
        _ (float): Unused.

    Returns:
        None
    """
    # Generate test data with no duplicates
    arr1, arr2 = call_func(size, subset_ratio, 0)

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

    Args:
        call_func (TestFuncCallable): The function that generates test data.
        size (int): The size of the arrays to generate.
        subset_ratio (float): The ratio of elements in `arr1` that are also in `arr2`.
        _ (float): Unused.

    Returns:
        None
    """

    arr1, arr2 = call_func(size, subset_ratio, 0)
    # Check if arr1 is a strict subset of arr2
    assert au.is_strict_subset(arr2, arr1)

    arr1, arr2 = call_func(size, 1, 0)
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

    Args:
        call_func (TestFuncCallable): The function that generates test data.
        size (int): The size of the arrays to generate.
        subset_ratio (float): Unused.
        _ (float): Unused.

    Returns:
        None
    """
    # Generate test data where arr1 is a strict superset of arr2
    arr1, arr2 = call_func(size, subset_ratio, 0)
    # Check if arr1 is a strict superset of arr2
    assert au.is_strict_superset(arr1, arr2)

    # Generate test data where arr1 is equal to arr2
    arr1, arr2 = call_func(size, 1, 0)
    # Check if arr1 is not a strict superset of arr2
    assert not au.is_strict_superset(arr1, arr2)
