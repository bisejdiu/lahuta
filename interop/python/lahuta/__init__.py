"""High-level Python package entry for Lahuta."""

import sys

if sys.version_info < (3, 10):
    raise RuntimeError(f"Lahuta requires Python >= 3.10 (found {sys.version.split()[0]}).")

_missing_dependencies: list[str] = []
try:
    import numpy  # noqa: F401
except ImportError:
    _missing_dependencies.append("numpy")

try:
    import orjson  # noqa: F401
except ImportError:
    _missing_dependencies.append("orjson")

if _missing_dependencies:
    deps_str = ", ".join(_missing_dependencies)
    print(
        f"ERROR: Missing required dependencies: {deps_str}. Install with: pip install {' '.join(_missing_dependencies)}",
        file=sys.stderr,
    )
    _cpp_bindings_available = False
    _import_error = f"Missing dependencies: {deps_str}"
else:
    try:
        # fmt: off
        from .lib import lahuta as lxx
        from .lib.lahuta import ArpeggioContactsEngine, AtomRec, AtomType, Category, Contact, \
            AtomTypingMethod, ContactProvider, ContactSet, EntityID, EntityResolver, FastNS, KDIndex, \
            FeatureGroup, Flavor, GroupRec, IR, IdentityAnalyzerLuni, InteractionType, Kind, \
            LahutaSystem, LahutaSystemProperties, Logger, LuniFileProcessor, LuniPropertyResult, \
            MolStarContactsEngine, GetContactsEngine, NSResults, PropertyAnalyzerLuni, PropertyKey, \
            PropertyQueryLuni, Residue, Residues, RingRec, SearchOptions, Topology, \
            TopologyBuildingOptions, TopologyComputers, compute_angles, factorize, find_contacts, process_files

        rdkit     = lxx.rdkit
        metrics   = lxx.metrics
        neighbors = lxx.neighbors

        # must be kept here after the above import to not mess up cold access times
        from .config import logging as logging
        from .neighbors import NearestNeighbors

        # So `import lahuta.logging` and `from lahuta.logging import LogLevel` work
        sys.modules.setdefault(__name__ + ".logging", logging)

        _cpp_bindings_available = True
        _import_error = None
    except ImportError as e:
        _cpp_bindings_available = False
        _import_error = str(e)

        print(
            f"ERROR: Lahuta could not be imported: {e}\nEnsure Lahuta has been built and installed correctly. See INSTALL.md for details.",
            file=sys.stderr,
        )


def verify_bindings():
    return _cpp_bindings_available


def get_import_error():
    return _import_error if not _cpp_bindings_available else None


# fmt: off
if _cpp_bindings_available:
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
        "verify_bindings", "get_import_error", "logging", "rdkit",
    ]
else:
    # Fallback if C++ bindings are not available
    # NOTE: We probably want to throw if we fail to import the bindings. There is no
    # reason to use the Python package without the C++ core.
    __all__ = ["verify_bindings", "get_import_error"]

def _quick_self_test() -> bool:
    """Internal preflight and warm-up of array path."""
    if not _cpp_bindings_available:
        return False

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
    if _cpp_bindings_available:
        # best-effort, ignore failure
        _quick_self_test()
except Exception:
    pass
