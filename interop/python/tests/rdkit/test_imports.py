import importlib


def test_rdkit_submodule_import():
    import lahuta  # noqa: F401

    # check the rdkit submodule can be imported directly
    m = importlib.import_module("lahuta.rdkit")
    assert hasattr(m, "Conformer")

    # direct symbol import works
    from lahuta.rdkit import Conformer  # noqa: F401
