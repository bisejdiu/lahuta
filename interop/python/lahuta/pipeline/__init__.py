# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     print(string.Template("$a$b@$c").substitute(a="besian", b="sejdiu", c="gmail.com"))
#
from __future__ import annotations

from lahuta.lib import lahuta as _lib

from .core import ErrorInfo, ErrorRecord, Pipeline, PipelineResult, PyTaskFn
from .params import SystemParams, TopologyParams
from .tasks import ContactTask
from .types import FileOutput, InMemoryPolicy, OutputFormat, PipelineContext, ShardedOutput

DataField = _lib.pipeline.DataField

StageManager = _lib.pipeline.StageManager
ReportingLevel = _lib.pipeline.ReportingLevel

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
    "StageManager",
    "ReportingLevel",
    "ContactTask",
    "SystemParams",
    "TopologyParams",
    "DataField",
)
