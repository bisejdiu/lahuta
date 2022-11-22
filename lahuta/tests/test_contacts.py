import warnings
from pathlib import Path

import numpy as np
import pytest

import lahuta.contacts.contacts as contacts
from lahuta.core.neighbors import NeighborPairs
from lahuta.core.universe import Universe


class ExpectedResults:
    COVALENT = (
        np.array(
            [[805, 1174], [1204, 238], [916, 880], [1196, 198], [1174, 252], [43, 26]]
        ),
        (6, 2),
        np.array(
            [2.36713798, 1.81648471, 2.0938709, 1.81355522, 1.96590266, 1.32155373]
        ),
    )
    METALIC = (
        np.array([[805, 1174], [1174, 252]]),
        (2, 2),
        np.array([2.36713798, 1.96590266]),
    )
    CARBONYL = (
        np.array(
            [[802, 751], [759, 801], [236, 290], [998, 1039], [129, 195], [666, 632]]
        ),
        (8, 2),
        np.array(
            [
                [
                    3.32255232,
                    3.40379212,
                    3.3972972,
                    3.54279681,
                    3.24019192,
                    3.59262434,
                    3.57549888,
                    3.55704924,
                ]
            ]
        ),
    )
    HBOND = (
        np.array(
            [[749, 700], [749, 689], [844, 716], [844, 700], [851, 819], [697, 632]]
        ),
        (64, 2),
        np.array(
            [2.91437828, 3.22367449, 2.89316307, 3.2123208, 2.85801188, 2.87266228]
        ),
    )
    WEAK_HBOND = (
        np.array(
            [[651, 982], [651, 701], [916, 30], [916, 114], [583, 530], [985, 702]]
        ),
        (44, 2),
        np.array(
            [3.50017918, 3.44197959, 3.84160327, 3.55128948, 3.62837438, 3.54509803]
        ),
    )
    IONIC = (
        np.array(
            [
                [721, 942],
                [1159, 447],
                [1189, 515],
                [1189, 514],
                [1189, 513],
                [1189, 512],
            ]
        ),
        (11, 2),
        np.array(
            [3.95734407, 3.99556483, 3.99903327, 3.12490787, 3.38146656, 3.70707066]
        ),
    )
    AROMATIC = (
        np.array(
            [
                [1199, 250],
                [1199, 252],
                [1198, 250],
                [1198, 252],
                [1183, 252],
                [1183, 251],
            ]
        ),
        (20, 2),
        np.array(
            [3.42738844, 3.60225579, 2.95220069, 2.7691934, 3.73579025, 3.67549289]
        ),
    )
    HYDROPHOBIC = (
        np.array(
            [[760, 805], [760, 803], [793, 1211], [717, 982], [635, 1184], [635, 1186]]
        ),
        (150, 2),
        np.array(
            [4.48367049, 3.78603908, 4.08985156, 4.0888844, 4.08697226, 4.36636773]
        ),
    )
    POLAR_HBOND = (
        np.array(
            [[740, 667], [749, 700], [749, 689], [664, 632], [844, 819], [844, 716]]
        ),
        (101, 2),
        np.array(
            [3.0323696, 2.91437828, 3.22367449, 2.90402245, 3.31475641, 2.89316307]
        ),
    )
    POLAR_WEAK_HBOND = (
        np.array(
            [[759, 803], [759, 800], [667, 739], [651, 701], [632, 703], [610, 1002]]
        ),
        (95, 2),
        np.array(
            [3.48957345, 3.33343581, 3.2203678, 3.44197959, 3.40424396, 3.43800592]
        ),
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
    assert cov.pairs.shape == ExpectedResults.COVALENT[1]
    assert np.all(cov.pairs == ExpectedResults.COVALENT[0])
    assert np.allclose(cov.distances, ExpectedResults.COVALENT[2])


def test_metalic_neighbors(neighbors):
    met = contacts.metalic_neighbors(neighbors)
    assert met.pairs.shape == ExpectedResults.METALIC[1]
    assert np.all(met.pairs == ExpectedResults.METALIC[0])
    assert np.allclose(met.distances, ExpectedResults.METALIC[2])


def test_carbonyl_neighbors(neighbors):
    carb = contacts.carbonyl_neighbors(neighbors)
    assert np.all(carb.pairs.shape == ExpectedResults.CARBONYL[1])
    assert np.all(carb.pairs[:6] == ExpectedResults.CARBONYL[0])
    assert np.allclose(carb.distances, ExpectedResults.CARBONYL[2])


def test_hbond_neighbors(neighbors):
    hb = contacts.hbond_neighbors(neighbors)
    assert np.all(hb.pairs.shape == ExpectedResults.HBOND[1])
    assert np.all(hb.pairs[:6] == ExpectedResults.HBOND[0])
    assert np.allclose(hb.distances[:6], ExpectedResults.HBOND[2])


def test_weak_hbond_neighbors(neighbors):
    whb = contacts.weak_hbond_neighbors(neighbors)
    assert np.all(whb.pairs.shape == ExpectedResults.WEAK_HBOND[1])
    assert np.all(whb.pairs[:6] == ExpectedResults.WEAK_HBOND[0])
    assert np.allclose(whb.distances[:6], ExpectedResults.WEAK_HBOND[2])


def test_ionic_neighbors(neighbors):
    ionic = contacts.ionic_neighbors(neighbors)
    assert np.all(ionic.pairs.shape == ExpectedResults.IONIC[1])
    assert np.all(ionic.pairs[:6] == ExpectedResults.IONIC[0])
    assert np.allclose(ionic.distances[:6], ExpectedResults.IONIC[2])


def test_aromatic_neighbors(neighbors):
    aromatic = contacts.aromatic_neighbors(neighbors)
    assert np.all(aromatic.pairs.shape == ExpectedResults.AROMATIC[1])
    assert np.all(aromatic.pairs[:6] == ExpectedResults.AROMATIC[0])
    assert np.allclose(aromatic.distances[:6], ExpectedResults.AROMATIC[2])


def test_hydrophobic_neighbors(neighbors):
    hydrophobic = contacts.hydrophobic_neighbors(neighbors)
    assert np.all(hydrophobic.pairs.shape == ExpectedResults.HYDROPHOBIC[1])
    assert np.all(hydrophobic.pairs[:6] == ExpectedResults.HYDROPHOBIC[0])
    assert np.allclose(hydrophobic.distances[:6], ExpectedResults.HYDROPHOBIC[2])


def test_polar_hbond_neighbors(neighbors):
    phb = contacts.polar_hbond_neighbors(neighbors)
    assert np.all(phb.pairs.shape == ExpectedResults.POLAR_HBOND[1])
    assert np.all(phb.pairs[:6] == ExpectedResults.POLAR_HBOND[0])
    assert np.allclose(phb.distances[:6], ExpectedResults.POLAR_HBOND[2])


def test_weak_polar_hbond_neighbors(neighbors):
    pwhb = contacts.weak_polar_hbond_neighbors(neighbors)
    assert np.all(pwhb.pairs.shape == ExpectedResults.POLAR_WEAK_HBOND[1])
    assert np.all(pwhb.pairs[:6] == ExpectedResults.POLAR_WEAK_HBOND[0])
    assert np.allclose(pwhb.distances[:6], ExpectedResults.POLAR_WEAK_HBOND[2])
