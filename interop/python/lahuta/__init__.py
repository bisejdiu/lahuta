"""High-level Python package entry for Lahuta."""

import sys
from importlib.util import find_spec
from importlib.metadata import version, PackageNotFoundError

if sys.version_info < (3, 10):
    raise RuntimeError(f"Lahuta requires Python >= 3.10 (found {sys.version.split()[0]}).")


def _require_importable(dist: str, *, min_version: str | None = None) -> None:
    if find_spec(dist) is None:
        hint = f"pip install '{dist}{'>=' + min_version if min_version else ''}'"
        raise ImportError(f"Lahuta requires {dist}. Install it with: {hint}")

    if min_version:
        try:
            installed = version(dist)
        except PackageNotFoundError as e:
            raise ImportError(f"Lahuta requires {dist} >= {min_version}.") from e
        from packaging.version import Version

        if Version(installed) < Version(min_version):
            raise ImportError(f"Lahuta requires {dist} >= {min_version} (found {installed}).")


_require_importable("numpy", min_version="2.2")
_require_importable("orjson", min_version="3.11")

_missing_dependencies: list[str] = []
if find_spec("cloudpickle") is None:
    _missing_dependencies.append("cloudpickle")

if _missing_dependencies:
    deps_str = ", ".join(_missing_dependencies)
    raise ImportError(f"Lahuta requires {deps_str}. Install with: pip install {' '.join(_missing_dependencies)}")

try:
    # fmt: off
    from .lib import lahuta as lxx
    from .lib.lahuta import ArpeggioContactsEngine, AtomRec, AtomType, Category, Contact, \
        AtomTypingMethod, ContactProvider, ContactSet, EntityID, EntityResolver, FastNS, KDIndex, \
        FeatureGroup, Flavor, GroupRec, IR, IdentityAnalyzerLuni, InteractionType, Kind, \
        LahutaSystem, LahutaSystemProperties, Logger, LuniFileProcessor, LuniPropertyResult, \
        MolStarContactsEngine, GetContactsEngine, NSResults, PropertyAnalyzerLuni, PropertyKey, \
        PropertyQueryLuni, Residue, Residues, RingRec, SearchOptions, Topology, \
        TopologyBuildingOptions, TopologyComputers, compute_angles, factorize, find_contacts, process_files, vdw_radius

    rdkit     = lxx.rdkit
    metrics   = lxx.metrics
    neighbors = lxx.neighbors

    # must be kept here after the above import to not mess up cold access times
    from .config import logging as logging
    from .neighbors import NearestNeighbors

    # So `import lahuta.logging` and `from lahuta.logging import LogLevel` work
    sys.modules.setdefault(__name__ + ".logging", logging)

except ImportError as e:
    raise ImportError(
        f"Lahuta C++ bindings could not be imported: {e}\n"
        "Ensure Lahuta has been built and installed correctly. See INSTALL.md for details."
    ) from e


# fmt: off
__all__ = [
    'ArpeggioContactsEngine', 'AtomRec', 'AtomType', 'Category', 'Contact',
    'AtomTypingMethod', 'ContactProvider', 'ContactSet',
    'EntityID', 'EntityResolver', 'FastNS', 'KDIndex', 'FeatureGroup',
    'Flavor', 'GroupRec', 'IR', 'IdentityAnalyzerLuni', 'InteractionType', 'Kind',
    'LahutaSystem', 'LahutaSystemProperties', 'Logger', 'LuniFileProcessor',
    'LuniPropertyResult', 'MolStarContactsEngine', 'GetContactsEngine', 'NSResults',
    'PropertyAnalyzerLuni', 'PropertyKey', 'PropertyQueryLuni',
    'Residue', 'Residues', 'RingRec', 'SearchOptions', 'Topology',
    'TopologyBuildingOptions', 'TopologyComputers', 'compute_angles',
    'factorize', 'find_contacts', 'metrics', 'neighbors', 'process_files', 'NearestNeighbors',
    "vdw_radius", "logging", "rdkit",
]

def _quick_self_test() -> bool:
    """Internal preflight and warm-up of array path."""
    try:
        from .lib.lahuta import LahutaSystem as ___LahutaSystem
        ___rdkit = rdkit
        ___mol = ___rdkit.RWMol()
        _ = ___mol.addAtom(6)
        ___conf = ___rdkit.Conformer(1)
        ___conf.set3D(True)
        ___conf.setAtomPos(0, ___rdkit.Point3D(0.0, 0.0, 0.0))
        ___mol.addConformer(___conf, True)
        __sys = ___LahutaSystem.create(___mol)
        arr = __sys.props.positions_view  # warm and validate
        ok = getattr(arr, "shape", None) == (1, 3)
        del __sys, ___mol, ___conf
        return bool(ok)
    except Exception:
        return False

#
# Import-time warm-up: prime numpy ndarray and pybind11/RDKit casters so all
# numpy-backed arrays (positions, neighbors, distances) avoid cold start cost.
# A microoptimization, if there ever was one.     - Besian, September 2025
#
try:
    # best-effort, ignore failure
    _quick_self_test()
except Exception:
    pass
