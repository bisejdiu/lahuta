import json
import warnings
from pathlib import Path

import numpy as np
import pytest

import lahuta.contacts as C
from lahuta.contacts import F
from lahuta.contacts.atom_plane import AtomPlaneContacts
from lahuta.contacts.plane_plane import PlanePlaneContacts
from lahuta.core.universe import Universe

# pylint: disable=redefined-outer-name


class ExpectedResults:
    with open(Path(__file__).parent / "data" / "1KX2.json", "r", encoding="utf-8") as f:
        data = json.load(f)

    COVALENT = data["COVALENT"]
    METALIC = data["METALIC"]
    CARBONYL = data["CARBONYL"]
    HBOND = data["HBOND"]
    WEAK_HBOND = data["WEAK_HBOND"]
    IONIC = data["IONIC"]
    AROMATIC = data["AROMATIC"]
    HYDROPHOBIC = data["HYDROPHOBIC"]
    POLAR_HBOND = data["POLAR_HBOND"]
    WEAK_POLAR_HBOND = data["WEAK_POLAR_HBOND"]
    VDW = data["VDW"]
    CARBONPI = data["CARBONPI"]
    CATIONPI = data["CATIONPI"]
    DONORPI = data["DONORPI"]
    SULPHURPI = data["SULPHURPI"]
    PLANEPLANE = data["PLANEPLANE"]


@pytest.fixture(scope="session")
def data_loader():
    """
    Fixture to load the data for the tests.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        path_obj = Path(__file__).parent / "data" / "1KX2.pdb"
        universe = Universe(str(path_obj))
    ns = universe.compute_neighbors()
    return ns


@pytest.fixture(scope="session")
def neighbors(data_loader):
    """Helper fixture to get neighbor pairs."""
    return data_loader


@pytest.fixture(scope="session")
def atom_plane(data_loader):
    """Helper fixture to get atomplane."""
    atomplane = AtomPlaneContacts(data_loader)
    return atomplane


@pytest.fixture(scope="session")
def plane_plane(data_loader):
    """Helper fixture to get planeplane."""
    planeplane = PlanePlaneContacts(data_loader)

    return planeplane


@pytest.mark.parametrize(
    "neighbor_func, expected_result",
    [
        (F.covalent_neighbors, ExpectedResults.COVALENT),
        (F.metalic_neighbors, ExpectedResults.METALIC),
        (F.carbonyl_neighbors, ExpectedResults.CARBONYL),
        (F.hbond_neighbors, ExpectedResults.HBOND),
        (F.weak_hbond_neighbors, ExpectedResults.WEAK_HBOND),
        (F.ionic_neighbors, ExpectedResults.IONIC),
        (F.aromatic_neighbors, ExpectedResults.AROMATIC),
        (F.hydrophobic_neighbors, ExpectedResults.HYDROPHOBIC),
        (F.polar_hbond_neighbors, ExpectedResults.POLAR_HBOND),
        (F.weak_polar_hbond_neighbors, ExpectedResults.WEAK_POLAR_HBOND),
        (F.vdw_neighbors, ExpectedResults.VDW),
        (F.sulphur_pi, ExpectedResults.SULPHURPI),
        (F.carbon_pi, ExpectedResults.CARBONPI),
        (F.cation_pi, ExpectedResults.CATIONPI),
        (F.donor_pi, ExpectedResults.DONORPI),
        (F.plane_plane_neighbors, ExpectedResults.PLANEPLANE),
    ],
)
def test_atom_atom_neighbor_funcs(neighbor_func, expected_result, neighbors):
    """Test the neighbors."""
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
    "contact_class, expected_result",
    [
        (C.CovalentContacts, ExpectedResults.COVALENT),
        (C.MetalicContacts, ExpectedResults.METALIC),
        (C.CarbonylContacts, ExpectedResults.CARBONYL),
        (C.HBondContacts, ExpectedResults.HBOND),
        (C.WeakHBondContacts, ExpectedResults.WEAK_HBOND),
        (C.IonicContacts, ExpectedResults.IONIC),
        (C.AromaticContacts, ExpectedResults.AROMATIC),
        (C.HydrophobicContacts, ExpectedResults.HYDROPHOBIC),
        (C.PolarHBondContacts, ExpectedResults.POLAR_HBOND),
        (C.WeakPolarHBondContacts, ExpectedResults.WEAK_POLAR_HBOND),
        (C.VanDerWaalsContacts, ExpectedResults.VDW),
        (C.SulphurPi, ExpectedResults.SULPHURPI),
        (C.CarbonPi, ExpectedResults.CARBONPI),
        (C.CationPi, ExpectedResults.CATIONPI),
        (C.DonorPi, ExpectedResults.DONORPI),
    ],
)
def test_atom_atom_neighbor_classes(contact_class, expected_result, neighbors):
    """Test the neighbors."""
    instance = contact_class(neighbors)
    # instance.run()
    result = instance.results

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
    "contact_func, expected_result",
    [
        (
            lambda ap: ap.cation_pi(),
            ExpectedResults.CATIONPI,
        ),
        (
            lambda ap: ap.donor_pi(),
            ExpectedResults.DONORPI,
        ),
        (
            lambda ap: ap.sulphur_pi(),
            ExpectedResults.SULPHURPI,
        ),
        (
            lambda ap: ap.carbon_pi(),
            ExpectedResults.CARBONPI,
        ),
    ],
)
def test_atomplane_contacts(contact_func, expected_result, atom_plane):
    """Test the contacts."""
    result = contact_func(atom_plane)
    pairs, distances = np.array(result.pairs[:6]), np.array(result.distances[:6])

    expected_pairs = np.array(expected_result["pairs"])

    # Reshape the expected_pairs to match the shape of the pairs
    if expected_pairs.size == 0:
        expected_pairs = expected_pairs.reshape((0, 2))

    assert result.pairs.shape[0] == expected_result["shapex"]
    assert np.all(pairs == expected_pairs)
    assert np.allclose(distances, expected_result["distances"], atol=1e-3)


def test_planeplane_neighbors(plane_plane):
    """Test the plane_plane neighbors."""
    plane_ns = plane_plane.results
    pairs, distances = np.array(plane_ns.pairs), np.array(plane_ns.distances)
    assert plane_ns.pairs.shape[0] == ExpectedResults.PLANEPLANE["shapex"]
    assert np.all(pairs == ExpectedResults.PLANEPLANE["pairs"])
    assert np.allclose(distances, ExpectedResults.PLANEPLANE["distances"], atol=1e-3)
