import gzip
import json
import warnings
from pathlib import Path
from typing import Callable, Type, Union

import numpy as np
import pytest
from _pytest.fixtures import FixtureRequest

import lahuta.contacts as C
from lahuta import Luni
from lahuta.contacts import F
from lahuta.contacts.atom_plane import AtomPlaneContacts
from lahuta.contacts.base import ContactAnalysis
from lahuta.contacts.plane_plane import PlanePlaneContacts
from lahuta.core.neighbors import NeighborPairs
from lahuta.tests import X2, DNABound, Rhodopsin


ContactFunction = Callable[[NeighborPairs], NeighborPairs]
ContactDict = dict[str, Union[list[int], int]]

pytestmark = pytest.mark.contacts

FILE_PAIRS = [("1kx2.json.gz", X2()), ("1gzm.json.gz", Rhodopsin()), ("3q2y.json.gz", DNABound())]


@pytest.fixture(scope="session", params=FILE_PAIRS)
def data_loader(request: FixtureRequest) -> tuple[NeighborPairs, dict[str, ContactDict]]:
    json_file, file_obj = request.param

    # Load ExpectedResults from the JSON file
    with gzip.open(Path(__file__).parent / "data" / "results" / json_file, "rt", encoding="utf-8") as file:
        data = json.load(file)

    # Load universe from the pdb_file
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        universe = Luni(str(file_obj))

    ns = universe.compute_neighbors(res_dif=1)
    return ns, data


@pytest.fixture(scope="session")
def neighbors(data_loader: tuple[NeighborPairs, dict[str, ContactDict]]) -> NeighborPairs:
    """Helper fixture to get neighbor pairs."""
    return data_loader[0]


@pytest.fixture(scope="session")
def atom_plane(data_loader: tuple[NeighborPairs, dict[str, ContactDict]]) -> AtomPlaneContacts:
    """Helper fixture to get atomplane."""
    atomplane = AtomPlaneContacts(data_loader[0])
    return atomplane


@pytest.fixture(scope="session")
def plane_plane(data_loader: tuple[NeighborPairs, dict[str, ContactDict]]) -> PlanePlaneContacts:
    """Helper fixture to get planeplane."""
    planeplane = PlanePlaneContacts(data_loader[0])

    return planeplane


@pytest.mark.parametrize(
    "neighbor_func, expected_key",
    [
        (F.covalent_neighbors, "COVALENT"),
        (F.metalic_neighbors, "METALIC"),
        (F.carbonyl_neighbors, "CARBONYL"),
        (F.hbond_neighbors, "HBOND"),
        (F.weak_hbond_neighbors, "WEAK_HBOND"),
        (F.ionic_neighbors, "IONIC"),
        (F.aromatic_neighbors, "AROMATIC"),
        (F.hydrophobic_neighbors, "HYDROPHOBIC"),
        (F.polar_hbond_neighbors, "POLAR_HBOND"),
        (F.weak_polar_hbond_neighbors, "WEAK_POLAR_HBOND"),
        (F.vdw_neighbors, "VDW"),
        (F.sulphur_pi, "SULPHURPI"),
        (F.carbon_pi, "CARBONPI"),
        (F.cation_pi, "CATIONPI"),
        (F.donor_pi, "DONORPI"),
        (F.plane_plane_neighbors, "PLANEPLANE"),
    ],
)
def test_atom_atom_neighbor_funcs(
    neighbor_func: Callable[[NeighborPairs], NeighborPairs],
    expected_key: str,
    neighbors: NeighborPairs,
    data_loader: tuple[NeighborPairs, dict[str, ContactDict]],
) -> None:
    """Test the neighbors."""

    _, expected_results = data_loader
    expected_result = expected_results[expected_key]

    result = neighbor_func(neighbors)
    pairs, distances = result.pairs, result.distances

    expected_pairs = np.array(expected_result["pairs"])

    # # Reshape the expected_pairs to match the shape of the pairs
    if expected_pairs.size == 0:
        expected_pairs = expected_pairs.reshape((0, 2))

    assert result.pairs.shape[1] == 2
    assert result.pairs.shape[0] == expected_result["shapex"]
    assert np.all(pairs == expected_pairs)
    assert np.allclose(distances.tolist(), expected_result["distances"], atol=1e-3)


@pytest.mark.parametrize(
    "contact_class, expected_key",
    [
        (C.CovalentContacts, "COVALENT"),
        (C.MetalicContacts, "METALIC"),
        (C.CarbonylContacts, "CARBONYL"),
        (C.HBondContacts, "HBOND"),
        (C.WeakHBondContacts, "WEAK_HBOND"),
        (C.IonicContacts, "IONIC"),
        (C.AromaticContacts, "AROMATIC"),
        (C.HydrophobicContacts, "HYDROPHOBIC"),
        (C.PolarHBondContacts, "POLAR_HBOND"),
        (C.WeakPolarHBondContacts, "WEAK_POLAR_HBOND"),
        (C.VanDerWaalsContacts, "VDW"),
        (C.SulphurPi, "SULPHURPI"),
        (C.CarbonPi, "CARBONPI"),
        (C.CationPi, "CATIONPI"),
        (C.DonorPi, "DONORPI"),
    ],
)
def test_atom_atom_neighbor_classes(
    contact_class: Type[ContactAnalysis],
    expected_key: str,
    neighbors: NeighborPairs,
    data_loader: tuple[NeighborPairs, dict[str, ContactDict]],
) -> None:
    """Test the neighbors."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        instance = contact_class(neighbors)
    # instance.run()
    result = instance.results

    _, expected_results = data_loader
    expected_result = expected_results[expected_key]

    pairs, distances = result.pairs, result.distances

    expected_pairs = np.array(expected_result["pairs"])

    # Reshape the expected_pairs to match the shape of the pairs
    if expected_pairs.size == 0:
        expected_pairs = expected_pairs.reshape((0, 2))

    assert result.pairs.shape[1] == 2
    assert result.pairs.shape[0] == expected_result["shapex"]
    assert np.all(pairs == expected_pairs)
    assert np.allclose(distances, expected_result["distances"], atol=1e-3)


@pytest.mark.parametrize(
    "contact_func_name, expected_key",
    [
        ("cation_pi", "CATIONPI"),
        ("donor_pi", "DONORPI"),
        ("sulphur_pi", "SULPHURPI"),
        ("carbon_pi", "CARBONPI"),
    ],
)
def test_atomplane_contacts(
    contact_func_name: str,
    expected_key: str,
    atom_plane: AtomPlaneContacts,
    data_loader: tuple[NeighborPairs, dict[str, ContactDict]],
) -> None:
    """Test the contacts."""

    contact_func = getattr(atom_plane, contact_func_name)
    result = contact_func()
    pairs, distances = result.pairs, result.distances

    _, expected_results = data_loader
    expected_result = expected_results[expected_key]

    expected_pairs = np.array(expected_result["pairs"])

    # Reshape the expected_pairs to match the shape of the pairs
    if expected_pairs.size == 0:
        expected_pairs = expected_pairs.reshape((0, 2))

    assert result.pairs.shape[0] == expected_result["shapex"]
    assert np.all(pairs == expected_pairs)
    assert np.allclose(distances, expected_result["distances"], atol=1e-3)


def test_planeplane_neighbors(
    plane_plane: PlanePlaneContacts, data_loader: tuple[NeighborPairs, dict[str, ContactDict]]
) -> None:
    """Test the plane_plane neighbors."""
    plane_ns = plane_plane.results
    pairs, distances = plane_ns.pairs, plane_ns.distances

    # Get the expected result from data_loader
    _, expected_results = data_loader
    expected_plane_plane = expected_results["PLANEPLANE"]

    assert plane_ns.pairs.shape[0] == expected_plane_plane["shapex"]
    assert np.all(pairs == expected_plane_plane["pairs"])
    assert np.allclose(distances, expected_plane_plane["distances"], atol=1e-3)
