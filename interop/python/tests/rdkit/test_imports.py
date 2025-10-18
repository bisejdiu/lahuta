"""Tests that verify the rdkit submodule can be imported and basic symbols are available."""

import importlib

import pytest


def test_rdkit_submodule_import():
    try:
        import lahuta  # noqa: F401
    except ImportError:
        pytest.skip("C++ bindings not available! rdkit test skipped")

    # check the rdkit submodule can be imported directly
    m = importlib.import_module("lahuta.rdkit")
    assert hasattr(m, "Conformer")

    # direct symbol import works
    from lahuta.rdkit import Conformer  # noqa: F401
