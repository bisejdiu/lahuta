from __future__ import annotations

from lahuta.lib import lahuta as _lib

from .core import ErrorInfo, ErrorRecord, Pipeline, PipelineResult, PyTaskFn
from .params import SystemParams, TopologyParams
from .tasks import ContactTask
from .types import FileOutput, InMemoryPolicy, OutputFormat, PipelineContext, ShardedOutput

BackpressureConfig = _lib.pipeline.BackpressureConfig
OnFull = _lib.pipeline.OnFull
get_default_backpressure_config = _lib.pipeline.get_default_backpressure_config
set_default_backpressure_config = _lib.pipeline.set_default_backpressure_config
set_default_max_queue_bytes = _lib.pipeline.set_default_max_queue_bytes

__all__ = (
    "PyTaskFn",
    "ErrorRecord",
    "ErrorInfo",
    "Pipeline",
    "PipelineResult",
    "PipelineContext",
    "InMemoryPolicy",
    "OutputFormat",
    "FileOutput",
    "ShardedOutput",
    "BackpressureConfig",
    "OnFull",
    "get_default_backpressure_config",
    "set_default_backpressure_config",
    "set_default_max_queue_bytes",
    "ContactTask",
    "SystemParams",
    "TopologyParams",
)
