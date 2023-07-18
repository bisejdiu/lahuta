import json
import warnings
from pathlib import Path

import numpy as np
import pytest

import lahuta.contacts.contacts as contacts
from lahuta.contacts.plane import AtomPlaneContacts, PlanePlaneContacts
from lahuta.core.neighbors import NeighborPairs
from lahuta.core.universe import Universe

# pylint: disable=redefined-outer-name


class ExpectedResults:
    with open(Path(__file__).parent / "data" / "1KX2.json") as f:
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
        u = Universe(str(path_obj))
    n = u.compute_neighbors()
    return u, n


@pytest.fixture(scope="session")
def neighbors(data_loader):
    """Helper fixture to get neighbor pairs."""
    return data_loader[1]


@pytest.fixture(scope="session")
def ap(data_loader):
    """Helper fixture to get ap."""
    ap = AtomPlaneContacts(data_loader[0])
    ap.compute_contacts()
    return ap


@pytest.fixture(scope="session")
def planeplane(data_loader):
    """Helper fixture to get planeplane."""
    planeplane = PlanePlaneContacts(data_loader[0])
    planeplane.compute_contacts()
    return planeplane


@pytest.mark.parametrize(
    "neighbor_func, expected_result",
    [
        (contacts.covalent_neighbors, ExpectedResults.COVALENT),
        (contacts.metalic_neighbors, ExpectedResults.METALIC),
        (contacts.carbonyl_neighbors, ExpectedResults.CARBONYL),
        (contacts.hbond_neighbors, ExpectedResults.HBOND),
        (contacts.weak_hbond_neighbors, ExpectedResults.WEAK_HBOND),
        (contacts.ionic_neighbors, ExpectedResults.IONIC),
        (contacts.aromatic_neighbors, ExpectedResults.AROMATIC),
        (contacts.hydrophobic_neighbors, ExpectedResults.HYDROPHOBIC),
        (contacts.polar_hbond_neighbors, ExpectedResults.POLAR_HBOND),
        (contacts.weak_polar_hbond_neighbors, ExpectedResults.WEAK_POLAR_HBOND),
        (contacts.vdw_neighbors, ExpectedResults.VDW),
    ],
)
def test_atom_atom_neighbors(neighbor_func, expected_result, neighbors):
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
    "contact_func, expected_result",
    [
        (
            lambda ap: ap.cation_pi.contacts(ap.neighbors, ap.angles),
            ExpectedResults.CATIONPI,
        ),
        (
            lambda ap: ap.donor_pi.contacts(ap.neighbors, ap.angles),
            ExpectedResults.DONORPI,
        ),
        (
            lambda ap: ap.sulphur_pi.contacts(ap.neighbors, ap.angles),
            ExpectedResults.SULPHURPI,
        ),
        (
            lambda ap: ap.carbon_pi.contacts(ap.neighbors, ap.angles),
            ExpectedResults.CARBONPI,
        ),
    ],
)
def test_atomplane_contacts(contact_func, expected_result, ap):
    """Test the contacts."""
    result = contact_func(ap)
    pairs, distances = np.array(result.pairs[:6]), np.array(result.distances[:6])

    expected_pairs = np.array(expected_result["pairs"])

    # Reshape the expected_pairs to match the shape of the pairs
    if expected_pairs.size == 0:
        expected_pairs = expected_pairs.reshape((0, 2))

    assert result.pairs.shape[0] == expected_result["shapex"]
    assert np.all(pairs == expected_pairs)
    assert np.allclose(distances, expected_result["distances"], atol=1e-3)


def test_planeplane_neighbors(planeplane):
    """Test the planeplane neighbors."""
    pairs, distances = np.array(planeplane.pairs), np.array(planeplane.distances)
    assert planeplane.pairs.shape[0] == ExpectedResults.PLANEPLANE["shapex"]
    assert np.all(pairs == ExpectedResults.PLANEPLANE["pairs"])
    assert np.allclose(distances, ExpectedResults.PLANEPLANE["distances"], atol=1e-3)
