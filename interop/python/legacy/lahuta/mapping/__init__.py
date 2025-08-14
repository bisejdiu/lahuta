"""Lahuta mapping submodule for sequence and structure alignment."""

try:
    from ..lib import mappingxx
    from ..lib.mappingxx import (
        # Enums
        AlignType,
        TMScoreThrMode,
        SeqType,
        
        # FoldSeek classes
        FoldSeekOps,
        PrefilterOptions,
        MatcherResult,
        Matcher,
        
        # Alignment classes
        LahutaAlignerBase,
        ProcessingConfig,
        LahutaProcessor,
        SeqData,
        AlignerResults,
        LahutaAligner,
        
        # Topological equivalency classes
        LahutaMapper,
        TopologyMapper,
        ContactEquivKey,
        EquivalencyConfig,
        TopologicalEquivalency,
    )
    
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
        "LahutaMapper",
        "TopologyMapper",
        "ContactEquivKey",
        "EquivalencyConfig",
        "TopologicalEquivalency",
    ]
    
except ImportError as e:
    # Handle case where C++ extension is not built
    import warnings
    warnings.warn(f"Could not import mapping C++ extension: {e}")
    __all__ = []
