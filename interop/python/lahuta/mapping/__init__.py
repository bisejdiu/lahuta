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
