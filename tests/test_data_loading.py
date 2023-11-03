import warnings
from MDAnalysis import Universe

import pytest

import numpy as np

from lahuta import Luni
from lahuta.tests.base import BaseFile
from lahuta.tests import X2, Rhodopsin, DNABound

EXPECTED_DATA = {
    'x2': {
        'ids': np.arange(0, 1249),
        'resids': np.concatenate([np.arange(-3, 0), np.arange(1, 79), np.array([90])]),
        'element_count': {'C': 400, 'FE': 1, 'H': 610, 'N': 106, 'O': 124, 'S': 8},
        'resname_count': {
            'ALA': 122,
            'ARG': 24,
            'ASN': 56,
            'ASP': 84,
            'CYS': 40,
            'GLN': 17,
            'GLU': 60,
            'GLY': 42,
            'HEC': 75,
            'HIS': 34,
            'ILE': 38,
            'LEU': 76,
            'LYS': 177,
            'MET': 68,
            'PHE': 20,
            'PRO': 56,
            'SER': 44,
            'THR': 70,
            'TRP': 24,
            'TYR': 42,
            'VAL': 80,
        },
        'most_common_names': {
            'C': '81',
            'O': '81',
            'N': '81',
            'CA': '81',
            'H': '76',
            'HA': '75',
            'CB': '75',
            'HB3': '63',
            'HB2': '63',
            'CG': '43',
        },
    },
    'rhodopsin': {
        'ids': np.arange(0, 5792),
        'resids': np.concatenate(
            [
                np.arange(0, 327),
                np.array([330, 331, 332]),
                np.array([1332, 1333, 1334]),
                np.arange(1340, 1350),
                np.arange(2001, 2021),
            ]
        ),
        'element_count': {'C': 3891, 'N': 817, 'O': 1030, 'S': 52, 'ZN': 2},
        'resname_count': {
            'ACE': 6,
            'ALA': 260,
            'ARG': 154,
            'ASN': 240,
            'ASP': 68,
            'BMA': 44,
            'C8E': 182,
            'CYS': 120,
            'GLN': 198,
            'GLU': 272,
            'GLY': 176,
            'HIS': 120,
            'HOH': 40,
            'ILE': 352,
            'LDA': 48,
            'LEU': 432,
            'LYS': 180,
            'MET': 256,
            'NAG': 112,
            'PEF': 44,
            'PHE': 682,
            'PLM': 68,
            'PRO': 252,
            'RET': 40,
            'SER': 144,
            'THR': 322,
            'TRP': 140,
            'TYR': 432,
            'VAL': 406,
            'ZN': 2,
        },
        'most_common_names': {
            'O': '698',
            'CA': '660',
            'C': '658',
            'N': '658',
            'CB': '616',
            'CG': '368',
            'CD1': '206',
            'CD2': '174',
            'CG2': '148',
            'CD': '126',
        },
    },
    'dnabound': {
        'ids': np.arange(0, 2740),
        'resids': np.concatenate(
            [
                np.arange(-12, -1),
                np.arange(1, 278),
                np.arange(285, 298),
                np.array([3968]),
            ]
        ),
        'element_count': {'C': 1686, 'N': 443, 'O': 581, 'P': 22, 'S': 8},
        'resname_count': {
            'ALA': 70,
            'ARG': 126,
            'ASN': 64,
            'ASP': 112,
            'CYS': 6,
            'DA': 63,
            'DC': 152,
            'DG': 173,
            'DT': 80,
            'ET': 24,
            'GLN': 104,
            'GLU': 219,
            'GLY': 32,
            'GOL': 12,
            'HIS': 40,
            'HOH': 15,
            'ILE': 208,
            'LEU': 205,
            'LYS': 183,
            'MET': 56,
            'PHE': 132,
            'PRO': 84,
            'SER': 120,
            'THR': 112,
            'TYR': 264,
            'VAL': 84,
        },
        'most_common_names': {
            'O': '291',
            'C': '276',
            'CA': '276',
            'N': '276',
            'CB': '268',
            'CG': '164',
            'CD1': '85',
            'CD': '72',
            'CD2': '63',
            'CG2': '54',
        },
    },
}


def read_luni(pdb_file: str) -> Luni:
    """Read a PDB file and return a Luni object."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Luni(pdb_file)


def read_luni_from_mda(pdb_file: str) -> Luni:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        univ = Universe(pdb_file)
        return Luni(univ.atoms)


def run_per_file_tests(key: str, luni: Luni) -> None:
    assert luni.arc is not None
    assert luni.arc.residues is not None

    # atom ids
    assert np.array_equal(luni.arc.atoms.ids, EXPECTED_DATA[key]["ids"])
    assert np.array_equal(np.unique(luni.arc.residues.resids), EXPECTED_DATA[key]["resids"])

    # element counts
    unique, counts = np.unique(luni.arc.atoms.elements, return_counts=True)
    assert dict(zip(unique, counts)) == EXPECTED_DATA[key]["element_count"]

    # resname counts
    unique, counts = np.unique(luni.arc.residues.resnames, return_counts=True)
    assert dict(zip(unique, counts)) == EXPECTED_DATA[key]["resname_count"]

    # most common atom names
    unique, counts = np.unique(luni.arc.atoms.names, return_counts=True)
    top10 = np.asarray((unique, counts)).T[np.argsort(counts)][::-1][:10]
    most_common_names = dict(zip(top10[:, 0], top10[:, 1]))
    assert most_common_names == EXPECTED_DATA[key]["most_common_names"]


pdb_classes = [(X2, "x2"), (Rhodopsin, "rhodopsin"), (DNABound, "dnabound")]
TEST_PARAMS = [(name, cls, flag) for cls, name in pdb_classes for flag in [True, False]]


@pytest.mark.parametrize("pdb_name, pdb_class, pdb_flag", TEST_PARAMS)
def test_read_pdb(pdb_name: str, pdb_class: type[BaseFile], pdb_flag: bool) -> None:
    luni = read_luni(str(pdb_class(pdb_flag)))
    run_per_file_tests(pdb_name, luni)


@pytest.mark.parametrize("pdb_name, pdb_class, pdb_flag", TEST_PARAMS)
def test_univ_init(pdb_name: str, pdb_class: type[BaseFile], pdb_flag: bool) -> None:
    if not pdb_flag:
        return None
    luni = read_luni_from_mda(str(pdb_class(pdb_flag)))
    run_per_file_tests(pdb_name, luni)
