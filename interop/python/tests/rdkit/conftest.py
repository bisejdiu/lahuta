import importlib

import pytest


@pytest.fixture(scope="module")
def rdkit():
    return importlib.import_module("lahuta.rdkit")
