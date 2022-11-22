import json
import warnings
from pathlib import Path

import numpy as np
import pytest

import lahuta.contacts.contacts as contacts
from lahuta.core.neighbors import NeighborPairs
from lahuta.core.universe import Universe


# we will load the data from the json file
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


class DataLoader:
    def __init__(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.u = Universe(Path(__file__).parent / "data" / "1KX2.pdb")
        self.n = self.u.compute_neighbors()

    def __call__(self):
        return self.u, self.n


@pytest.fixture
def neighbors():
    return DataLoader().n


def test_covalent_neighbors(neighbors):
    cov = contacts.covalent_neighbors(neighbors)
    pairs, distances = np.array(cov.pairs[:6]), np.array(cov.distances[:6])
    assert cov.pairs.shape[0] == ExpectedResults.COVALENT["shapex"]
    assert np.all(pairs == ExpectedResults.COVALENT["pairs"])
    assert np.allclose(distances, ExpectedResults.COVALENT["distances"], atol=1e-3)


def test_metalic_neighbors(neighbors):
    met = contacts.metalic_neighbors(neighbors)
    pairs, distances = np.array(met.pairs), np.array(met.distances)
    assert met.pairs.shape[0] == ExpectedResults.METALIC["shapex"]
    assert np.all(np.array(met.pairs) == ExpectedResults.METALIC["pairs"])
    assert np.allclose(distances, ExpectedResults.METALIC["distances"], atol=1e-3)


def test_carbonyl_neighbors(neighbors):
    carb = contacts.carbonyl_neighbors(neighbors)
    pairs, distances = np.array(carb.pairs[:6]), np.array(carb.distances[:6])
    assert carb.pairs.shape[0] == ExpectedResults.CARBONYL["shapex"]
    assert np.all(pairs == ExpectedResults.CARBONYL["pairs"])
    assert np.allclose(distances, ExpectedResults.CARBONYL["distances"], atol=1e-3)


def test_hbond_neighbors(neighbors):
    hb = contacts.hbond_neighbors(neighbors)
    pairs, distances = np.array(hb.pairs[:6]), np.array(hb.distances[:6])
    assert hb.pairs.shape[0] == ExpectedResults.HBOND["shapex"]
    assert np.all(pairs == ExpectedResults.HBOND["pairs"])
    assert np.allclose(distances, ExpectedResults.HBOND["distances"], atol=1e-3)


def test_weak_hbond_neighbors(neighbors):
    whb = contacts.weak_hbond_neighbors(neighbors)
    pairs, distances = np.array(whb.pairs[:6]), np.array(whb.distances[:6])
    assert whb.pairs.shape[0] == ExpectedResults.WEAK_HBOND["shapex"]
    assert np.all(pairs == ExpectedResults.WEAK_HBOND["pairs"])
    assert np.allclose(distances, ExpectedResults.WEAK_HBOND["distances"], atol=1e-3)


def test_ionic_neighbors(neighbors):
    ionic = contacts.ionic_neighbors(neighbors)
    pairs, distances = np.array(ionic.pairs[:6]), np.array(ionic.distances[:6])
    assert ionic.pairs.shape[0] == ExpectedResults.IONIC["shapex"]
    assert np.all(pairs == ExpectedResults.IONIC["pairs"])
    assert np.allclose(distances, ExpectedResults.IONIC["distances"], atol=1e-3)


def test_aromatic_neighbors(neighbors):
    aromatic = contacts.aromatic_neighbors(neighbors)
    pairs, distances = np.array(aromatic.pairs[:6]), np.array(aromatic.distances[:6])
    assert aromatic.pairs.shape[0] == ExpectedResults.AROMATIC["shapex"]
    assert np.all(pairs == ExpectedResults.AROMATIC["pairs"])
    assert np.allclose(distances, ExpectedResults.AROMATIC["distances"])


def test_hydrophobic_neighbors(neighbors):
    hydrophobic = contacts.hydrophobic_neighbors(neighbors)
    pairs, distances = np.array(hydrophobic.pairs[:6]), np.array(
        hydrophobic.distances[:6]
    )
    assert hydrophobic.pairs.shape[0] == ExpectedResults.HYDROPHOBIC["shapex"]
    assert np.all(pairs == ExpectedResults.HYDROPHOBIC["pairs"])
    assert np.allclose(distances, ExpectedResults.HYDROPHOBIC["distances"])


def test_polar_hbond_neighbors(neighbors):
    phb = contacts.polar_hbond_neighbors(neighbors)
    pairs, distances = np.array(phb.pairs[:6]), np.array(phb.distances[:6])
    assert phb.pairs.shape[0] == ExpectedResults.POLAR_HBOND["shapex"]
    assert np.all(pairs == ExpectedResults.POLAR_HBOND["pairs"])
    assert np.allclose(distances, ExpectedResults.POLAR_HBOND["distances"])


def test_weak_polar_hbond_neighbors(neighbors):
    pwhb = contacts.weak_polar_hbond_neighbors(neighbors)
    pairs, distances = np.array(pwhb.pairs[:6]), np.array(pwhb.distances[:6])
    assert pwhb.pairs.shape[0] == ExpectedResults.WEAK_POLAR_HBOND["shapex"]
    assert np.all(np.array(pwhb.pairs[:6]) == ExpectedResults.WEAK_POLAR_HBOND["pairs"])
    assert np.allclose(distances, ExpectedResults.WEAK_POLAR_HBOND["distances"])


def test_vdw_neighbors(neighbors):
    vdw = contacts.vdw_neighbors(neighbors)
    pairs, distances = np.array(vdw.pairs[:6]), np.array(vdw.distances[:6])
    assert vdw.pairs.shape[0] == ExpectedResults.VDW["shapex"]
    assert np.all(pairs == ExpectedResults.VDW["pairs"])
    assert np.allclose(distances, ExpectedResults.VDW["distances"])
