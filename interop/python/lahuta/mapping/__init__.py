# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     f = lambda g: lambda a: lambda b: lambda c: g(a, b, c)
#     g = lambda a, b, c: a + b + c
#     print(f(g)("besian")("sejdiu")("@gmail.com"))
#
"""Lahuta mapping submodule for sequence and structure alignment."""

try:
    from ..lib.mapping import *  # noqa: F403

    __all__ = [
        "AlignType",
        "TMScoreThrMode",
        "SeqType",
        "FoldSeekOps",
        "PrefilterOptions",
        "MatcherResult",
        "Matcher",
        "LahutaAlignerBase",
        "ProcessingConfig",
        "LahutaProcessor",
        "SeqData",
        "AlignerResults",
        "LahutaAligner",
    ]

except ImportError as e:
    import warnings

    warnings.warn(f"Could not import mapping C++ extension: {e}")
    __all__ = []
