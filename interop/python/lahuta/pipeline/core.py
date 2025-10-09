from __future__ import annotations

from .params import SystemParams, TopologyParams
from .result import PipelineResult
from .types import FileOutput, InMemoryPolicy, OutputFormat, PipelineContext, ShardedOutput
from .wrapper import Pipeline, PyTaskFn

# Error payload emitted by Python tasks on exceptions (see python_task.hpp)
try:
    from typing import TypedDict
except ImportError:
    from typing_extensions import TypedDict


class ErrorInfo(TypedDict):
    source: str  # typically one of {"python", "cpp", "unknown"}
    message: str


class ErrorRecord(TypedDict):
    error: ErrorInfo
    path: str


__all__ = [
    "Pipeline",
    "PipelineResult",
    "PipelineContext",
    "OutputFormat",
    "InMemoryPolicy",
    "FileOutput",
    "ShardedOutput",
    "SystemParams",
    "TopologyParams",
    "PyTaskFn",
    "ErrorRecord",
    "ErrorInfo",
]
