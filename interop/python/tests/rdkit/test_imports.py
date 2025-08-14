"""Tests that verify the rdkit submodule can be imported and basic symbols are available."""

import importlib

import pytest


def test_rdkit_submodule_import():
    import lahuta

    if not lahuta.verify_bindings():
        pytest.skip("C++ bindings not available! rdkit test skipped")

    # check the rdkit submodule can be imported directly
    m = importlib.import_module("lahuta.rdkit")
    assert hasattr(m, "Conformer")

    # direct symbol import works
    from lahuta.rdkit import Conformer  # noqa: F401
