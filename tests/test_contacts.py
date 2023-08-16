import json
import warnings
from pathlib import Path
from typing import Callable, Dict, List, Tuple, Type, Union

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

# pylint: disable=redefined-outer-name
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

ContactFunction = Callable[[NeighborPairs], NeighborPairs]
ContactDict = Dict[str, Union[List[int], int]]

pytestmark = pytest.mark.contacts

file_pairs = [("1kx2.json", "1kx2.pdb"), ("1gzm.json", "1gzm.cif")]


@pytest.fixture(scope="session", params=file_pairs)
def data_loader(request: FixtureRequest) -> Tuple[NeighborPairs, Dict[str, ContactDict]]:
    json_file, pdb_file = request.param

    # Load ExpectedResults from the JSON file
    with open(Path(__file__).parent / "data" / "results" / json_file, "r", encoding="utf-8") as file:
        data = json.load(file)

    # Load universe from the pdb_file
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        path_obj = Path(__file__).parent / "data" / "results" / pdb_file
        universe = Luni(str(path_obj))

    ns = universe.compute_neighbors()
    return ns, data


@pytest.fixture(scope="session")
def neighbors(data_loader: Tuple[NeighborPairs, Dict[str, ContactDict]]) -> NeighborPairs:
    """Helper fixture to get neighbor pairs."""
    return data_loader[0]


@pytest.fixture(scope="session")
def atom_plane(data_loader: Tuple[NeighborPairs, Dict[str, ContactDict]]) -> AtomPlaneContacts:
    """Helper fixture to get atomplane."""
    atomplane = AtomPlaneContacts(data_loader[0])
    return atomplane


@pytest.fixture(scope="session")
def plane_plane(data_loader: Tuple[NeighborPairs, Dict[str, ContactDict]]) -> PlanePlaneContacts:
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
    data_loader: Tuple[NeighborPairs, Dict[str, ContactDict]],
) -> None:
    """Test the neighbors."""

    _, expected_results = data_loader
    expected_result = expected_results[expected_key]

    result = neighbor_func(neighbors)
    pairs, distances = np.array(result.pairs[:6]), np.array(result.distances[:6])

    expected_pairs = np.array(expected_result["pairs"])

    # Reshape the expected_pairs to match the shape of the pairs
    if expected_pairs.size == 0:
        expected_pairs = expected_pairs.reshape((0, 2))

    assert result.pairs.shape[1] == 2
    assert result.pairs.shape[0] == expected_result["shapex"]
    assert np.all(pairs == expected_pairs)
    assert np.allclose(distances, expected_result["distances"], atol=1e-3)


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
    data_loader: Tuple[NeighborPairs, Dict[str, ContactDict]],
) -> None:
    """Test the neighbors."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        instance = contact_class(neighbors)
    # instance.run()
    result = instance.results

    _, expected_results = data_loader
    expected_result = expected_results[expected_key]

    pairs, distances = np.array(result.pairs[:6]), np.array(result.distances[:6])

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
    data_loader: Tuple[NeighborPairs, Dict[str, ContactDict]],
) -> None:
    """Test the contacts."""

    contact_func = getattr(atom_plane, contact_func_name)
    result = contact_func()
    pairs, distances = np.array(result.pairs[:6]), np.array(result.distances[:6])

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
    plane_plane: PlanePlaneContacts, data_loader: Tuple[NeighborPairs, Dict[str, ContactDict]]
) -> None:
    """Test the plane_plane neighbors."""
    plane_ns = plane_plane.results
    pairs, distances = np.array(plane_ns.pairs), np.array(plane_ns.distances)

    # Get the expected result from data_loader
    _, expected_results = data_loader
    expected_plane_plane = expected_results["PLANEPLANE"]

    assert plane_ns.pairs.shape[0] == expected_plane_plane["shapex"]
    assert np.all(pairs[:6] == expected_plane_plane["pairs"])
    assert np.allclose(distances[:6], expected_plane_plane["distances"], atol=1e-3)
