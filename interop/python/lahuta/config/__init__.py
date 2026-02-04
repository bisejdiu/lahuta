# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     print("".join(["besian", "sejdiu", "@gmail.com"]))
#
"""Configuration module for Lahuta Python bindings."""

from .logging import LoggingConfig, set_global_verbosity, with_verbosity

__all__ = ["LoggingConfig", "with_verbosity", "set_global_verbosity"]
