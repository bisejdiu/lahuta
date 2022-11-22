import warnings
from pathlib import Path

import numpy as np
import pytest

import lahuta.contacts.contacts as contacts
from lahuta.core.neighbors import NeighborPairs
from lahuta.core.universe import Universe


class ExpectedResults:
    COVALENT = np.array(
        [[805, 1174], [1204, 238], [916, 880], [1196, 198], [1174, 252], [43, 26]]
    )
    METALIC = np.array([[805, 1174], [1174, 252]])
    CARBONYL = np.array(
        [[802, 751], [759, 801], [236, 290], [998, 1039], [129, 195], [666, 632]]
    )
    HBOND = np.array(
        [[749, 700], [749, 689], [844, 716], [844, 700], [851, 819], [697, 632]]
    )
    WEAK_HBOND = np.array(
        [[651, 982], [651, 701], [916, 30], [916, 114], [583, 530], [985, 702]]
    )
    IONIC = np.array(
        [[721, 942], [1159, 447], [1189, 515], [1189, 514], [1189, 513], [1189, 512]]
    )
    AROMATIC = np.array(
        [[1199, 250], [1199, 252], [1198, 250], [1198, 252], [1183, 252], [1183, 251]]
    )
    HYDROPHOBIC = np.array(
        [[760, 805], [760, 803], [793, 1211], [717, 982], [635, 1184], [635, 1186]]
    )
    POLAR_HBOND = np.array(
        [[740, 667], [749, 700], [749, 689], [664, 632], [844, 819], [844, 716]]
    )
    POLAR_WEAK_HBOND = np.array(
        [[759, 803], [759, 800], [667, 739], [651, 701], [632, 703], [610, 1002]]
    )


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
    assert np.all(cov.pairs == ExpectedResults.COVALENT)


def test_metalic_neighbors(neighbors):
    met = contacts.metalic_neighbors(neighbors)
    assert np.all(met.pairs == ExpectedResults.METALIC)


def test_carbonyl_neighbors(neighbors):
    carb = contacts.carbonyl_neighbors(neighbors)
    assert np.all(carb.pairs[:6] == ExpectedResults.CARBONYL)


def test_hbond_neighbors(neighbors):
    hb = contacts.hbond_neighbors(neighbors)
    assert np.all(hb.pairs[:6] == ExpectedResults.HBOND)


def test_weak_hbond_neighbors(neighbors):
    whb = contacts.weak_hbond_neighbors(neighbors)
    assert np.all(whb.pairs[:6] == ExpectedResults.WEAK_HBOND)


def test_ionic_neighbors(neighbors):
    ionic = contacts.ionic_neighbors(neighbors)
    assert np.all(ionic.pairs[:6] == ExpectedResults.IONIC)


def test_aromatic_neighbors(neighbors):
    aromatic = contacts.aromatic_neighbors(neighbors)
    assert np.all(aromatic.pairs[:6] == ExpectedResults.AROMATIC)


def test_hydrophobic_neighbors(neighbors):
    hydrophobic = contacts.hydrophobic_neighbors(neighbors)
    assert np.all(hydrophobic.pairs[:6] == ExpectedResults.HYDROPHOBIC)


def test_polar_hbond_neighbors(neighbors):
    phb = contacts.polar_hbond_neighbors(neighbors)
    assert np.all(phb.pairs[:6] == ExpectedResults.POLAR_HBOND)


def test_weak_polar_hbond_neighbors(neighbors):
    pwhb = contacts.weak_polar_hbond_neighbors(neighbors)
    assert np.all(pwhb.pairs[:6] == ExpectedResults.POLAR_WEAK_HBOND)
