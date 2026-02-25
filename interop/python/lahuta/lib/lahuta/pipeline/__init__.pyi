"""
Pipeline bindings
"""
from __future__ import annotations
import datetime
import lahuta.lib.lahuta
import typing
from . import sources
__all__: list[str] = ['BackpressureConfig', 'Block', 'DataField', 'DropLatest', 'DropOldest', 'LmdbSink', 'MemorySink', 'ModelPackTask', 'ModelPayload', 'NdjsonSink', 'OnFull', 'PipelineContext', 'ReportingLevel', 'ShardedNdjsonSink', 'Sink', 'SlowSink', 'StageManager', 'Task', 'get_default_backpressure_config', 'set_default_backpressure_config', 'set_default_max_queue_bytes', 'sources']
class BackpressureConfig:
    """
    Sink backpressure configuration
    """
    max_batch_bytes: int
    max_batch_msgs: int
    max_queue_bytes: int
    max_queue_msgs: int
    offer_wait_slice: datetime.timedelta
    on_full: OnFull
    required: bool
    writer_threads: int
    def __init__(self) -> None:
        ...
    def validate(self) -> None:
        ...
class DataField:
    Dssp: typing.ClassVar[DataField]  # value = <DataField.Dssp: 16>
    DsspView: typing.ClassVar[DataField]  # value = <DataField.DsspView: 256>
    Metadata: typing.ClassVar[DataField]  # value = <DataField.Metadata: 1>
    Plddt: typing.ClassVar[DataField]  # value = <DataField.Plddt: 8>
    PlddtView: typing.ClassVar[DataField]  # value = <DataField.PlddtView: 128>
    Positions: typing.ClassVar[DataField]  # value = <DataField.Positions: 4>
    PositionsView: typing.ClassVar[DataField]  # value = <DataField.PositionsView: 64>
    Sequence: typing.ClassVar[DataField]  # value = <DataField.Sequence: 2>
    SequenceView: typing.ClassVar[DataField]  # value = <DataField.SequenceView: 32>
    __members__: typing.ClassVar[dict[str, DataField]]  # value = {'Metadata': <DataField.Metadata: 1>, 'Sequence': <DataField.Sequence: 2>, 'Positions': <DataField.Positions: 4>, 'Plddt': <DataField.Plddt: 8>, 'Dssp': <DataField.Dssp: 16>, 'SequenceView': <DataField.SequenceView: 32>, 'PositionsView': <DataField.PositionsView: 64>, 'PlddtView': <DataField.PlddtView: 128>, 'DsspView': <DataField.DsspView: 256>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class LmdbSink(Sink):
    def __init__(self, db: ..., batch_size: int = 1024) -> None:
        ...
class MemorySink(Sink):
    def __init__(self) -> None:
        ...
    def clear(self) -> None:
        """
        Clear collected payloads.
        """
    def result(self) -> list[str]:
        """
        Return collected payloads.
        """
    def result_bytes(self) -> list:
        """
        Return collected payloads as bytes.
        """
class ModelPackTask(Task):
    def __init__(self, channel: str = 'db') -> None:
        ...
class ModelPayload:
    @property
    def dssp(self) -> typing.Any:
        ...
    @property
    def dssp_view(self) -> typing.Any:
        ...
    @property
    def metadata(self) -> typing.Any:
        ...
    @property
    def plddts(self) -> typing.Any:
        ...
    @property
    def plddts_view(self) -> typing.Any:
        ...
    @property
    def positions(self) -> typing.Any:
        ...
    @property
    def positions_view(self) -> typing.Any:
        ...
    @property
    def sequence(self) -> typing.Any:
        ...
    @property
    def sequence_view(self) -> typing.Any:
        ...
class NdjsonSink(Sink):
    def __init__(self, path: str) -> None:
        ...
    def file(self) -> str:
        ...
class OnFull:
    """
    Policy when a sink queue is full
    """
    Block: typing.ClassVar[OnFull]  # value = <OnFull.Block: 0>
    DropLatest: typing.ClassVar[OnFull]  # value = <OnFull.DropLatest: 1>
    DropOldest: typing.ClassVar[OnFull]  # value = <OnFull.DropOldest: 2>
    __members__: typing.ClassVar[dict[str, OnFull]]  # value = {'Block': <OnFull.Block: 0>, 'DropLatest': <OnFull.DropLatest: 1>, 'DropOldest': <OnFull.DropOldest: 2>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class PipelineContext:
    def frame_metadata(self) -> typing.Any:
        """
        Current frame metadata.
        """
    def get(self, arg0: str) -> typing.Any:
        ...
    def get_bytes(self, arg0: str) -> typing.Any:
        ...
    def get_json(self, arg0: str) -> typing.Any:
        ...
    def get_system(self) -> typing.Any:
        ...
    def get_text(self, arg0: str) -> typing.Any:
        ...
    def get_topology(self) -> typing.Any:
        ...
    def set_bytes(self, arg0: str, arg1: typing.Any) -> None:
        ...
    def set_json(self, arg0: str, arg1: typing.Any) -> None:
        ...
    def set_text(self, arg0: str, arg1: str) -> None:
        ...
    @property
    def conformer_id(self) -> int:
        """
        Return the conformer/frame identifier for the current item.
        """
    @property
    def model_payload(self) -> typing.Any:
        """
        Structured view over lazily loaded LMDB payload slices.
        """
    @property
    def path(self) -> str:
        """
        Return the input item path (e.g., file path).
        """
    @property
    def session_id(self) -> typing.Any:
        """
        Return the session identifier associated with this item.
        """
    @property
    def timestamp_ps(self) -> typing.Any:
        """
        Optional simulation timestamp in picoseconds.
        """
class ReportingLevel:
    BASIC: typing.ClassVar[ReportingLevel]  # value = <ReportingLevel.BASIC: 1>
    DEBUG: typing.ClassVar[ReportingLevel]  # value = <ReportingLevel.DEBUG: 2>
    OFF: typing.ClassVar[ReportingLevel]  # value = <ReportingLevel.OFF: 0>
    __members__: typing.ClassVar[dict[str, ReportingLevel]]  # value = {'OFF': <ReportingLevel.OFF: 0>, 'BASIC': <ReportingLevel.BASIC: 1>, 'DEBUG': <ReportingLevel.DEBUG: 2>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class ShardedNdjsonSink(Sink):
    @typing.overload
    def __init__(self, out_dir: str, shard_size: int) -> None:
        ...
    @typing.overload
    def __init__(self, out_dir: str, shard_size: int, max_shard_bytes: int) -> None:
        ...
    def files(self) -> list[str]:
        ...
class Sink:
    pass
class SlowSink(Sink):
    def __init__(self, sleep_seconds: float, flush_sleep_seconds: float = 0.0) -> None:
        """
        Test sink that sleeps during write()/flush. Used for flush-timeout testing.
        """
    def count(self) -> int:
        """
        Return number of writes completed.
        """
class StageManager:
    def __init__(self, source: sources.Source) -> None:
        ...
    def add_contacts(self, name: str, provider: lahuta.lib.lahuta.ContactProvider, interaction_types: typing.Any = None, channel: str | None = None, out_fmt: str = 'json', thread_safe: bool = True) -> None:
        ...
    def add_python(self, name: str, depends: list[str] = [], fn: typing.Any, channel: str | None = None, serialize: bool = True, store: bool = True) -> None:
        ...
    def add_task(self, name: str, depends: list[str] = [], task: Task, thread_safe: bool = True) -> None:
        ...
    def compile(self) -> None:
        ...
    def connect_sink(self, channel: str, sink: Sink, backpressure: BackpressureConfig | None = None) -> None:
        ...
    def describe_graph(self) -> list[tuple[str, list[str]]]:
        ...
    def get_auto_builtins(self) -> bool:
        ...
    def get_flush_timeout(self) -> float:
        ...
    def get_reporting_level(self) -> ReportingLevel:
        ...
    def get_system_params(self) -> dict:
        ...
    def get_task_data_fields(self, name: str) -> list[DataField]:
        ...
    def get_topology_params(self) -> dict:
        ...
    def last_run_report(self) -> typing.Any:
        ...
    def run(self, threads: int = 4, flush_timeout: float | None = None) -> None:
        ...
    def set_auto_builtins(self, arg0: bool) -> None:
        ...
    def set_flush_timeout(self, timeout: float) -> None:
        ...
    def set_reporting_level(self, arg0: ReportingLevel) -> None:
        ...
    def set_system_params(self, arg0: dict) -> None:
        ...
    def set_task_data_fields(self, name: str, fields: list[DataField]) -> None:
        ...
    def set_topology_params(self, arg0: dict) -> None:
        ...
    def sorted_tasks(self) -> list[str]:
        ...
    def stats(self) -> list:
        ...
class Task:
    pass
def get_default_backpressure_config() -> BackpressureConfig:
    """
    Return a copy of the default sink backpressure configuration.
    """
def set_default_backpressure_config(arg0: BackpressureConfig) -> None:
    """
    Set the default sink backpressure configuration used for new sinks.
    """
def set_default_max_queue_bytes(bytes: int) -> None:
    """
    Update only the byte capacity of the default sink queue (clamps batch bytes accordingly).
    """
Block: OnFull  # value = <OnFull.Block: 0>
DropLatest: OnFull  # value = <OnFull.DropLatest: 1>
DropOldest: OnFull  # value = <OnFull.DropOldest: 2>
