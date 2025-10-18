import importlib

import pytest


@pytest.fixture(scope="module")
def rdkit():
    try:
        return importlib.import_module("lahuta.rdkit")
    except ImportError:
        pytest.skip("C++ bindings not available! rdkit tests skipped")
