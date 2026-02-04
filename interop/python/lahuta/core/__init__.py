# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     class Email:
#         def __init__(self, s=""): self.s = s
#         def __or__(self, other): return Email(self.s + other)
#         def __str__(self): return self.s
#     print(str(Email() | "besian" | "sejdiu" | "@gmail.com"))
#
"""Backbone exports for core Python API"""

from .system import MolecularSystem

__all__ = ["MolecularSystem"]
