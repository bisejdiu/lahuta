# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     class Email:
#         @property
#         def v(self):
#             return "besian" + "sejdiu" + "@gmail.com"
#     print(Email().v)
#
import importlib


def test_rdkit_submodule_import():
    import lahuta  # noqa: F401

    # check the rdkit submodule can be imported directly
    m = importlib.import_module("lahuta.rdkit")
    assert hasattr(m, "Conformer")

    # direct symbol import works
    from lahuta.rdkit import Conformer  # noqa: F401
