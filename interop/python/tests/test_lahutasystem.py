"""
This test requires exact result generation. It verifies:
- Shapes/dtypes of numpy dtype properties
- Exact SHA-256 digests for all primary arrays
- Topology build pattern & accessors
- Neighbor-search counts and summary stats
"""

from __future__ import annotations

import hashlib

import numpy as np
import pytest

import lahuta as lxx


# fmt: off
def sha256_of_ndarray(arr: np.ndarray, name: str = "") -> str:
    """SHA256 over array metadata + content."""
    if not isinstance(arr, np.ndarray):
        raise TypeError(f"{name or 'array'} is not a NumPy array")

    h = hashlib.sha256()
    header = f"LHX1|{name}|{arr.dtype.str}|{arr.shape}".encode("utf-8")
    h.update(header)

    if arr.dtype.kind in ("U", "S") or arr.dtype == object:
        flat = arr.ravel(order="C").tolist()
        sep = "\x1f" 
        payload = sep.join("" if x is None else str(x) for x in flat).encode("utf-8")
        h.update(payload)
    else:
        c = np.ascontiguousarray(arr)
        h.update(c.view(np.uint8).tobytes())

    return h.hexdigest()


EXPECTED = {
    "n_atoms_file": 671,
    "n_atoms_filtered": 344,
    "centroid": np.array([-4.21122355, 4.16459314, -1.54197914], dtype=np.float64),
    # Hashes for arrays (shape, dtype, sha256)
    "hashes": {
        "indices":     ((671,), "int32",  "1ab38760f81915e948aa53672172de2fec931eeabaf0f3efeb66e2cbf4313d21"),
        "atom_nums":   ((671,), "int32",  "111a012c1b720b50c803036c1b4b9cd9fb4c7e276b6e1e0034da795e0551fc3d"),
        "resids":      ((671,), "int32",  "2a50c31b68221ed954d0442a60a09e6e29135d552bd442a5f08f76f5d33c2159"),
        "resindices":  ((671,), "int32",  "695ecbe032a769c5bfab939faf8f08e54c294cf6e1ed204637c1c41d8de0fede"),
        "names":       ((671,), "object", "a6a3a9d28124a3bbbb57b64e025703bdf19566bc0ee2982692c42611df8f0ecb"),
        "symbols":     ((671,), "object", "07468f9a8f5a0a4869a92a87a985b57393dd9bf89e58ae113779a553e3201905"),
        "elements":    ((671,), "object", "6ef6cc55261792e42927e0e76bbbf2492e3eeaf723fe276b7bc7d74c2ba0e008"),
        "resnames":    ((671,), "object", "d9b3997afaf24028ab05b62d7eb21a7702ed0389fc12dc4942364c458eab51b0"),
        "chainlabels": ((671,), "object", "f083ba8173908dd36daf76e2739d7f13848ca7fa28189d59acb455a8bfff3473"),
        "positions":   ((671, 3), "float64", "fffa7095f5dd892e4e9be8ecfe2bfe4e91d22fe9b0ffdcb691a9fea1adbd751a"),
    },

    "neighbors_file": {
        "cutoff"  : 4.5,
        "res_dif" : 1,
        "count"   : 1500,
        "mean"    : 3.957,
        "min"     : 2.584,
        "max"     : 4.498,
        "rtol"    : 5e-4,
        "atol"    : 1e-6,
    },

    "neighbors_filtered": {
        "cutoff"  : 4.0,
        "res_dif" : 0,
        "count"   : 1706,
        "mean"    : 2.747,
        "min"     : 1.224,
        "max"     : 3.999,
        "rtol"    : 5e-4,
        "atol"    : 1e-6,
    },
}


@pytest.fixture(scope="session")
def props(luni: lxx.LahutaSystem) -> lxx.LahutaSystemProperties:
    return luni.props


@pytest.fixture(scope="session")
def filtered_backbone(luni: lxx.LahutaSystem) -> lxx.LahutaSystem:
    names = luni.props.names
    keep = [i for i, nm in enumerate(names) if str(nm).strip() in {"N", "CA", "C", "O"}]
    return luni.filter(keep)

def test_basic_shapes_dtypes_and_hashes(luni: lxx.LahutaSystem, props: lxx.LahutaSystemProperties) -> None:
    assert luni.n_atoms == EXPECTED["n_atoms_file"]

    arrays: dict[str, np.ndarray] = {
        "indices":     props.indices,
        "atom_nums":   props.atom_nums,
        "resids":      props.resids,
        "resindices":  props.resindices,
        "names":       props.names,
        "symbols":     props.symbols,
        "elements":    props.elements,
        "resnames":    props.resnames,
        "chainlabels": props.chainlabels,
        "positions":   props.positions,
    }

    expected = EXPECTED["hashes"]
    for key, arr in arrays.items():
        exp_shape, exp_dtype, exp_hash = expected[key]
        assert tuple(arr.shape) == exp_shape, f"{key}: shape mismatch"
        assert str(arr.dtype) == exp_dtype, f"{key}: dtype mismatch (got {arr.dtype})"
        got_hash = sha256_of_ndarray(arr, name=key)
        assert got_hash == exp_hash, f"{key}: sha256 mismatch"


def test_centroid(props: lxx.LahutaSystemProperties) -> None:
    centroid = props.positions.mean(axis=0)
    assert centroid.shape == (3,)
    assert np.allclose(centroid, EXPECTED["centroid"], rtol=5e-6, atol=1e-8)


def test_topology_build_and_accessors(luni: lxx.LahutaSystem) -> None:
    # Configure and build
    luni.enable_only(lxx.TopologyComputers.Standard)
    luni.set_search_cutoff_for_bonds(1.9)
    assert luni.build_topology() is True
    assert luni.has_topology_built() is True

    # Execute Rings
    if not luni.is_computation_enabled(lxx.TopologyComputers.Rings):
        luni.enable_computation(lxx.TopologyComputers.Rings, True)
    assert luni.execute_computation(lxx.TopologyComputers.Rings) is True

    # Access RDKit handles
    mol  = luni.get_molecule()
    conf = luni.get_conformer()
    xyz  = luni.props.positions
    assert mol.getNumAtoms()  == luni.n_atoms
    assert conf.getNumAtoms() == luni.n_atoms
    assert xyz.shape == (luni.n_atoms, 3)

    # Topology handle should be retrievable
    topo = luni.get_topology()
    assert topo is not None


def _sqrt_distances(ns) -> np.ndarray:
    if getattr(ns, "get_sqrt_distances", None) is not None:
        return ns.get_sqrt_distances()
    dv = getattr(ns, "distances_view", None)
    if dv is not None:
        return np.sqrt(dv)
    return np.sqrt(ns.distances)


def test_neighbor_search_file_system(luni: lxx.LahutaSystem) -> None:
    cfg = EXPECTED["neighbors_file"]
    ns = luni.find_neighbors(cutoff=cfg["cutoff"], res_dif=cfg["res_dif"])
    ns = ns.filter(cfg["cutoff"])
    dij = _sqrt_distances(ns)
    assert dij.shape[0] == cfg["count"]

    mean = float(dij.mean())
    mn   = float(dij.min())
    mx   = float(dij.max())
    assert np.isclose(mean, cfg["mean"], rtol=cfg["rtol"], atol=cfg["atol"])
    assert np.isclose(mn,   cfg["min"],  rtol=cfg["rtol"], atol=cfg["atol"])
    assert np.isclose(mx,   cfg["max"],  rtol=cfg["rtol"], atol=cfg["atol"])


def test_neighbor_search_filtered_system(filtered_backbone: lxx.LahutaSystem) -> None:
    cfg = EXPECTED["neighbors_filtered"]
    ns  = filtered_backbone.find_neighbors(cutoff=cfg["cutoff"], res_dif=cfg["res_dif"])
    ns  = ns.filter(cfg["cutoff"])
    dij = _sqrt_distances(ns)
    assert filtered_backbone.n_atoms == EXPECTED["n_atoms_filtered"]
    assert dij.shape[0] == cfg["count"]

    mean = float(dij.mean())
    mn   = float(dij.min())
    mx   = float(dij.max())
    assert np.isclose(mean, cfg["mean"], rtol=cfg["rtol"], atol=cfg["atol"])
    assert np.isclose(mn,   cfg["min"],  rtol=cfg["rtol"], atol=cfg["atol"])
    assert np.isclose(mx,   cfg["max"],  rtol=cfg["rtol"], atol=cfg["atol"])


def test_hash_stability_under_basic_copies(props: lxx.LahutaSystemProperties) -> None:
    hyp = pytest.importorskip("hypothesis",            reason="Hypothesis not installed")
    st  = pytest.importorskip("hypothesis.strategies", reason="Hypothesis not installed")

    arrays = {
        "indices"     : props.indices,
        "atom_nums"   : props.atom_nums,
        "resids"      : props.resids,
        "resindices"  : props.resindices,
        "names"       : props.names,
        "symbols"     : props.symbols,
        "elements"    : props.elements,
        "resnames"    : props.resnames,
        "chainlabels" : props.chainlabels,
        "positions"   : props.positions,
    }
    keys = list(arrays.keys())

    @hyp.given(st.sampled_from(keys), st.sampled_from(["copy", "contig", "tolist_roundtrip"]))
    def _prop(key: str, mode: str) -> None:
        arr = arrays[key]
        h0 = sha256_of_ndarray(arr, name=key)

        if mode == "copy":
            arr2 = arr.copy()
        elif mode == "contig":
            arr2 = np.ascontiguousarray(arr)
        elif mode == "tolist_roundtrip":
            # Re-materialize the same data & dtype. For object/string this keeps exact str() values
            arr2 = np.array(arr.tolist(), dtype=arr.dtype, copy=True)
        else:
            raise AssertionError("unreachable")

        h1 = sha256_of_ndarray(arr2, name=key)
        assert h1 == h0

    _prop()
