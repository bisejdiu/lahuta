"""Tests for basic Luni functionality."""
import hashlib

import numpy as np
import pytest
from numpy.typing import NDArray

from lahuta import Luni
from lahuta.utils.download_files import X2, BaseFile

TEST_PARAMS = [(X2, "x2")]


def read_luni(pdb_file: str) -> Luni:
    """Read a PDB file and return a Luni object."""
    return Luni(pdb_file)


@pytest.fixture(scope="module")  # or use scope="session"
def luni() -> Luni:
    """Return a Luni object."""
    pdb_class = X2
    luni = read_luni(pdb_class().file_loc)
    return luni


def unique_size_counts(arr: NDArray[np.str_ | np.int32]) -> tuple[NDArray[np.str_ | np.int32], NDArray[np.int32]]:
    """Compute the unique values and their counts in an array."""
    unique, counts = np.unique(arr, return_counts=True)
    return unique, counts


def unique_size_sum_counts_sum(arr: NDArray[np.str_ | np.int32]) -> tuple[int, int]:
    """Compute the unique values and their counts in an array."""
    unique, counts = np.unique(arr, return_counts=True)
    return unique.sum(), counts.sum()


def unique_join_counts_sum(arr: NDArray[np.str_ | np.int32]) -> tuple[str, int]:
    """Compute the unique values and their counts in an array."""
    unique, counts = np.unique(arr, return_counts=True)
    unique_joined = "".join(unique)
    return unique_joined, counts.sum()


def hash_array(arr: NDArray[np.int32 | np.float32]) -> str:
    """Hash an array."""
    return hashlib.md5(arr.view(np.uint8).tobytes()).hexdigest()


@pytest.mark.parametrize(("pdb_class", "pdb_name"), TEST_PARAMS)
def test_attributes(luni: Luni, pdb_class: type[BaseFile], pdb_name: str) -> None:  # noqa: ARG001
    """Test the attributes of the Luni object."""
    assert luni.n_atoms == 1249
    assert luni.n_residues == 82
    assert luni.n_chains == 1
    assert luni.arc is not None

    chainauths, chainids, chainids = luni.chainauths, luni.chainids, luni.chainids
    assert chainauths is not None and chainids is not None and chainids is not None
    assert unique_size_counts(chainauths) == (np.array(["A"], dtype="<U10"), np.array([1249]))
    assert unique_size_counts(chainids) == (np.array([1]), np.array([1249]))
    assert unique_size_counts(chainids) == (np.array([1]), np.array([1249]))

    resids, resnames = luni.resids, luni.resnames
    assert resids is not None and resnames is not None
    assert unique_size_sum_counts_sum(resids) == (3165, 1249)
    assert unique_join_counts_sum(resnames) == ("ALAARGASNASPCYSGLNGLUGLYHECHISILELEULYSMETPHEPROSERTHRTRPTYRVAL", 1249)

    names, elements, types = luni.names, luni.elements, luni.types
    assert names is not None and elements is not None and types is not None
    assert unique_join_counts_sum(names) == (
        (
            "CC1AC1BC1CC1DC2AC2BC2CC2DC3AC3BC3CC3DC4AC4BC4CC4DCACAACABCACCADCBCBACBBCBCCBDCDCD1"
            "CD2CECE1CE2CE3CGCG1CG2CGACGDCH2CHACHBCHCCHDCMACMBCMCCMDCZCZ2CZ3FEHH1H2H3HAHA2HA3HA"
            "A1HAA2HABHACHAD1HAD2HBHB1HB2HB3HBA1HBA2HBB1HBB2HBB3HBC1HBC2HBC3HBD1HBD2HD1HD11HD12"
            "HD13HD2HD21HD22HD23HD3HEHE1HE2HE21HE22HE3HGHG1HG11HG12HG13HG2HG21HG22HG23HG3HHHH11"
            "HH12HH2HH21HH22HHAHHBHHCHHDHMA1HMA2HMA3HMB1HMB2HMB3HMC1HMC2HMC3HMD1HMD2HMD3HZHZ1HZ"
            "2HZ3NNANBNCNDND1ND2NENE1NE2NH1NH2NZOO1AO1DO2AO2DOD1OD2OE1OE2OGOG1OHOXTSDSG"
        ),
        1249,
    )
    assert unique_join_counts_sum(elements) == ("CFEHNOS", 1249)
    assert unique_join_counts_sum(types) == ("CFEHNOS", 1249)


def test_filtering() -> None:
    """Test the filtering of the Luni object."""
    luni = read_luni(X2().file_loc)

    assert all(luni.remove_ions().indices == luni.indices)  # X2 has no ions
    no_ligs = np.setdiff1d(luni.indices, luni.remove_ligands().indices)
    no_prot = np.setdiff1d(luni.indices, luni.remove_protein().indices)

    assert hash_array(no_ligs) == "253b0bb218a8ee56c1e4977c1c8740c3"
    assert hash_array(no_prot) == "9d58638694df48737464746d7d50ac0e"
    assert hash_array(luni.coordinates) == "27b1cafaf849b38ee70c9a8fae3303e3"


def test_luni_copy() -> None:
    """Test the copy of the Luni object."""
    luni = read_luni(X2().file_loc)

    mda, mda_copy = luni.to("mda"), luni.filter("protein").copy().to("mda")
    assert mda.n_atoms == mda_copy.universe.atoms.n_atoms == 1249
    assert mda.n_residues == mda_copy.universe.atoms.n_residues == 82
    assert mda.n_segments == mda_copy.universe.atoms.n_segments == 1

    mda, mda_deepcopy = luni.to("mda"), luni.filter("protein").deepcopy().to("mda")
    assert (mda.n_atoms == 1249) and (mda_deepcopy.n_atoms == 1174)
    assert (mda.n_residues == 82) and (mda_deepcopy.n_residues == 81)
    # TODO (bisejdiu): This gives an error
    # assert (mda.n_segments == 1) and (mda_deepcopy.n_segments == 1)


def test_extend_topology() -> None:
    """Test the extension of the topology of the Luni object."""
    luni = read_luni(X2().file_loc)
    rng = np.random.default_rng()
    rand_values = rng.random(luni.n_atoms)
    luni.extend_topology("topology_extend_test", rand_values)

    assert np.array_equal(luni.to("mda").topology_extend_test, rand_values)  # type: ignore[attr-defined]
