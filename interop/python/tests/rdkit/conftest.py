import importlib

import pytest


@pytest.fixture(scope="module")
def rdkit():
    import lahuta

    if not lahuta.verify_bindings():
        pytest.skip("C++ bindings not available! rdkit tests skipped")
    return importlib.import_module("lahuta.rdkit")
