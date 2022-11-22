import json
import warnings
from pathlib import Path

import numpy as np
import pytest

import lahuta.contacts.contacts as contacts
from lahuta.contacts.plane import AtomPlaneContacts, PlanePlaneContacts
from lahuta.core.neighbors import NeighborPairs
from lahuta.core.universe import Universe


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


class DataLoader:
    """
    Class to load the data for the tests.
    """

    def __init__(self):
        """Load the data for the tests."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.u = Universe(Path(__file__).parent / "data" / "1KX2.pdb")
        self.n = self.u.compute_neighbors()

    def __call__(self):
        """Return the data."""
        return self.u, self.n


@pytest.fixture
def neighbors():
    """Helper fixture to get neighbor pairs."""
    return DataLoader().n


@pytest.fixture
def atomplane():
    """Helper fixture to get atomplane."""
    atomplane = AtomPlaneContacts(DataLoader().u)
    atomplane.compute_contacts()
    return atomplane


@pytest.fixture
def planeplane():
    """Helper fixture to get planeplane."""
    planeplane = PlanePlaneContacts(DataLoader().u)
    planeplane.compute_contacts()
    return planeplane


def test_covalent_neighbors(neighbors):
    """Test the covalent neighbors."""
    cov = contacts.covalent_neighbors(neighbors)
    pairs, distances = np.array(cov.pairs[:6]), np.array(cov.distances[:6])
    assert cov.pairs.shape[0] == ExpectedResults.COVALENT["shapex"]
    assert np.all(pairs == ExpectedResults.COVALENT["pairs"])
    assert np.allclose(distances, ExpectedResults.COVALENT["distances"], atol=1e-3)


def test_metalic_neighbors(neighbors):
    """Test the metalic neighbors."""
    met = contacts.metalic_neighbors(neighbors)
    pairs, distances = np.array(met.pairs), np.array(met.distances)
    assert met.pairs.shape[0] == ExpectedResults.METALIC["shapex"]
    assert np.all(np.array(met.pairs) == ExpectedResults.METALIC["pairs"])
    assert np.allclose(distances, ExpectedResults.METALIC["distances"], atol=1e-3)


def test_carbonyl_neighbors(neighbors):
    """Test the carbonyl neighbors."""
    cbnyl = contacts.carbonyl_neighbors(neighbors)
    pairs, distances = np.array(cbnyl.pairs[:6]), np.array(cbnyl.distances[:6])
    assert cbnyl.pairs.shape[0] == ExpectedResults.CARBONYL["shapex"]
    assert np.all(pairs == ExpectedResults.CARBONYL["pairs"])
    assert np.allclose(distances, ExpectedResults.CARBONYL["distances"], atol=1e-3)


def test_ionic_neighbors(neighbors):
    """Test the ionic neighbors."""
    ionic = contacts.ionic_neighbors(neighbors)
    pairs, distances = np.array(ionic.pairs[:6]), np.array(ionic.distances[:6])
    assert ionic.pairs.shape[0] == ExpectedResults.IONIC["shapex"]
    assert np.all(pairs == ExpectedResults.IONIC["pairs"])
    assert np.allclose(distances, ExpectedResults.IONIC["distances"], atol=1e-3)


def test_aromatic_neighbors(neighbors):
    """Test the aromatic neighbors."""
    aromatic = contacts.aromatic_neighbors(neighbors)
    pairs, distances = np.array(aromatic.pairs[:6]), np.array(aromatic.distances[:6])
    assert aromatic.pairs.shape[0] == ExpectedResults.AROMATIC["shapex"]
    assert np.all(pairs == ExpectedResults.AROMATIC["pairs"])
    assert np.allclose(distances, ExpectedResults.AROMATIC["distances"])


def test_hydrophobic_neighbors(neighbors):
    """Test the hydrophobic neighbors."""
    hydrophobic = contacts.hydrophobic_neighbors(neighbors)
    pairs, distances = np.array(hydrophobic.pairs[:6]), np.array(
        hydrophobic.distances[:6]
    )
    assert hydrophobic.pairs.shape[0] == ExpectedResults.HYDROPHOBIC["shapex"]
    assert np.all(pairs == ExpectedResults.HYDROPHOBIC["pairs"])
    assert np.allclose(distances, ExpectedResults.HYDROPHOBIC["distances"])


def test_vdw_neighbors(neighbors):
    """Test the vdw neighbors."""
    vdw = contacts.vdw_neighbors(neighbors)
    pairs, distances = np.array(vdw.pairs[:6]), np.array(vdw.distances[:6])
    assert vdw.pairs.shape[0] == ExpectedResults.VDW["shapex"]
    assert np.all(pairs == ExpectedResults.VDW["pairs"])
    assert np.allclose(distances, ExpectedResults.VDW["distances"])


def test_hbond_neighbors(neighbors):
    """Test the hbond neighbors."""
    hb = contacts.hbond_neighbors(neighbors)
    pairs, distances = np.array(hb.pairs[:6]), np.array(hb.distances[:6])
    assert hb.pairs.shape[0] == ExpectedResults.HBOND["shapex"]
    assert np.all(pairs == ExpectedResults.HBOND["pairs"])
    assert np.allclose(distances, ExpectedResults.HBOND["distances"], atol=1e-3)


def test_weak_hbond_neighbors(neighbors):
    """Test the weak hbond neighbors."""
    whb = contacts.weak_hbond_neighbors(neighbors)
    pairs, distances = np.array(whb.pairs[:6]), np.array(whb.distances[:6])
    assert whb.pairs.shape[0] == ExpectedResults.WEAK_HBOND["shapex"]
    assert np.all(pairs == ExpectedResults.WEAK_HBOND["pairs"])
    assert np.allclose(distances, ExpectedResults.WEAK_HBOND["distances"], atol=1e-3)


def test_polar_hbond_neighbors(neighbors):
    """Test the polar hbond neighbors."""
    phb = contacts.polar_hbond_neighbors(neighbors)
    pairs, distances = np.array(phb.pairs[:6]), np.array(phb.distances[:6])
    assert phb.pairs.shape[0] == ExpectedResults.POLAR_HBOND["shapex"]
    assert np.all(pairs == ExpectedResults.POLAR_HBOND["pairs"])
    assert np.allclose(distances, ExpectedResults.POLAR_HBOND["distances"])


def test_weak_polar_hbond_neighbors(neighbors):
    """Test the weak polar hbond neighbors."""
    pwhb = contacts.weak_polar_hbond_neighbors(neighbors)
    pairs, distances = np.array(pwhb.pairs[:6]), np.array(pwhb.distances[:6])
    assert pwhb.pairs.shape[0] == ExpectedResults.WEAK_POLAR_HBOND["shapex"]
    assert np.all(pairs == ExpectedResults.WEAK_POLAR_HBOND["pairs"])
    assert np.allclose(distances, ExpectedResults.WEAK_POLAR_HBOND["distances"])


def test_carbonpi_neighbors(atomplane):
    """Test the carbonpi neighbors."""
    cpi = atomplane.carbon_pi.contacts(atomplane.neighbors, atomplane.angles)
    pairs, distances = np.array(cpi.pairs[:6]), np.array(cpi.distances[:6])
    assert cpi.pairs.shape[0] == ExpectedResults.CARBONPI["shapex"]
    assert np.all(pairs == ExpectedResults.CARBONPI["pairs"])
    assert np.allclose(distances, ExpectedResults.CARBONPI["distances"], atol=1e-3)


def test_cationpi_neighbors(atomplane):
    """Test the cationpi neighbors."""
    cpi = atomplane.cation_pi.contacts(atomplane.neighbors, atomplane.angles)
    assert cpi.pairs.shape[0] == ExpectedResults.CATIONPI["shapex"]
    assert np.all([] == ExpectedResults.CATIONPI["pairs"])
    assert np.allclose([], ExpectedResults.CATIONPI["distances"], atol=1e-3)


def test_donorpi_neighbors(atomplane):
    """Test the donorpi neighbors."""
    dpi = atomplane.donor_pi.contacts(atomplane.neighbors, atomplane.angles)
    pairs, distances = np.array(dpi.pairs[:6]), np.array(dpi.distances[:6])
    assert dpi.pairs.shape[0] == ExpectedResults.DONORPI["shapex"]
    assert np.all(pairs == ExpectedResults.DONORPI["pairs"])
    assert np.allclose(distances, ExpectedResults.DONORPI["distances"], atol=1e-3)


def test_sulphurpi_neighbors(atomplane):
    """Test the sulphurpi neighbors."""
    spi = atomplane.sulphur_pi.contacts(atomplane.neighbors, atomplane.angles)
    pairs, distances = np.array(spi.pairs[:6]), np.array(spi.distances[:6])
    assert spi.pairs.shape[0] == ExpectedResults.SULPHURPI["shapex"]
    assert np.all(pairs == ExpectedResults.SULPHURPI["pairs"])
    assert np.allclose(distances, ExpectedResults.SULPHURPI["distances"], atol=1e-3)


def test_planeplane_neighbors(planeplane):
    """Test the planeplane neighbors."""
    pairs, distances = np.array(planeplane.pairs), np.array(planeplane.distances)
    assert planeplane.pairs.shape[0] == ExpectedResults.PLANEPLANE["shapex"]
    assert np.all(pairs == ExpectedResults.PLANEPLANE["pairs"])
    assert np.allclose(distances, ExpectedResults.PLANEPLANE["distances"], atol=1e-3)
