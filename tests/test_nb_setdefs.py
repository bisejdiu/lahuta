from pathlib import Path
from typing import Tuple

import pytest
from _pytest.fixtures import FixtureRequest

import lahuta.utils.set_defs as sd
from lahuta import Universe
from lahuta.core.labeled_neighbors import LabeledNeighborPairs
from lahuta.msa.msa import MSAParser

# pylint: disable=redefined-outer-name
# pylint: disable=invalid-name
# pylint: disable=missing-function-docstring

pytestmark = pytest.mark.nb

# Type variables
T = Tuple[LabeledNeighborPairs, LabeledNeighborPairs]

file_triplets = [
    ('data/b2.cif', 'data/s5.cif', 'data/alig_b2s5.fasta'), 
]

# Fixtures
@pytest.fixture(params=file_triplets, ids=[f"triplet_{i}" for i in range(len(file_triplets))])
def setup_data(request: FixtureRequest) -> T:
    b2_data, s5_data, fasta_data = request.param
    b2_data_path = Path(__file__).parent / b2_data
    s5_data_path = Path(__file__).parent / s5_data
    fasta_data_path = Path(__file__).parent / fasta_data

    b2u = Universe(str(b2_data_path))
    s5u = Universe(str(s5_data_path))
    ns_b2 = b2u.compute_neighbors(res_dif=2)
    ns_s5 = s5u.compute_neighbors(res_dif=2)
    parser = MSAParser(str(fasta_data_path))
    seq_id_b2 = parser.get_seq_ids()[0]
    seq_b2 = parser.sequences[seq_id_b2]
    seq_id_s5 = parser.get_seq_ids()[1]
    seq_s5 = parser.sequences[seq_id_s5]

    s1 = ns_b2.map(seq_b2)
    s2 = ns_s5.map(seq_s5)

    return s1, s2

@pytest.fixture
def non_transformed_maps(setup_data: T) -> T:
    s1, s2 = setup_data
    return s1, s2

params = [
    (['ASP'], ['ASP', 'LEU'], 'atom_names'),
    (['ASP'], ['ASP', 'LEU'], 'resnames'),
    (['ALA', 'LEU', 'ILE', 'VAL'], ['ALA', 'LEU', 'ILE', 'VAL'], 'atom_names'),
]

@pytest.fixture(params=params)
def select_maps(setup_data: T, request: FixtureRequest) -> T:
    s1, s2 = setup_data
    transform_args = request.param
    res1, res2, label = transform_args
    return s1.select(resnames=res1).inverse().remove(label), s2.select(resnames=res2).remove(label)

@pytest.fixture(params=params)
def exclude_maps(setup_data: T, request: FixtureRequest) -> T:
    s1, s2 = setup_data
    transform_args = request.param
    res1, res2, label = transform_args
    return s1.exclude(resnames=res1).inverse().remove(label), s2.exclude(resnames=res2).remove(label)

@pytest.fixture(params=params)
def select_exclude_maps(setup_data: T, request: FixtureRequest) -> T:
    s1, s2 = setup_data
    transform_args = request.param
    res1, res2, label = transform_args
    return s1.select(resnames=res1).inverse().remove(label), s2.exclude(resnames=res2).remove(label)


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

# Test functions
def test_ops(non_transformed_maps: T) -> None:
    s1, s2 = non_transformed_maps
    run_common_tests(s1, s2)

def test_transformed_ops(select_maps: T) -> None:
    s1, s2 = select_maps
    run_common_tests(s1, s2)

def test_exclude_ops(exclude_maps: T) -> None:
    s1, s2 = exclude_maps
    run_common_tests(s1, s2)
