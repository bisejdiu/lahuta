# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     a, *rest = ["besian", "sejdiu", "@gmail.com"]
#     print(a + "".join(rest))
#
from __future__ import annotations

import time
from numbers import Integral
from pathlib import Path
from typing import (
    Any,
    Callable,
    Iterable,
    Literal,
    Protocol,
    Sequence,
    Type,
    TypeAlias,
    TypeVar,
    cast,
    overload,
)

from lahuta import logging
from lahuta.lib import lahuta as _lib
from lahuta.sources import (
    DatabaseHandleSource,
    DatabaseSource,
    LmdbSource,
    MdTrajectoriesSource,
    PipelineSource,
)

from ._discovery import _ast_infer_builtins_for_callable, _discover_builtins_for_callable
from ._inference import infer_python_task_format
from .params import SystemParams, TopologyParams
from .result import PipelineResult, _ChannelState
from .tasks import ContactTask
from .types import FileOutput, InMemoryPolicy, OutputFormat, PipelineContext, ShardedOutput

# fmt: off
DataField      = _lib.pipeline.DataField
StageManager   = _lib.pipeline.StageManager
ReportingLevel = _lib.pipeline.ReportingLevel
CppTask: TypeAlias = _lib.pipeline.Task

_DATABASE_SOURCE_TYPES: tuple[type[Any], ...] = tuple(
    t for t in (DatabaseSource, DatabaseHandleSource, LmdbSource) if t is not None
)

# Callable protocol for Python tasks that accept a PipelineContext and return any value
# that can be serialized into a single ndjson record. The return type is expressed
# by parameterizing PyTaskFn (e.g., PyTaskFn[dict[str, int]] or PyTaskFn[str]).
R_co    = TypeVar("R_co", covariant=True)
SchemaT = TypeVar("SchemaT")

class PyTaskFn(Protocol[R_co]):
    def __call__(self, __ctx: PipelineContext) -> R_co: ...


class Pipeline:
    """Pipeline builder and executor.

    A Pipeline orchestrates computational tasks. Tasks emit results to channels,
    which can be consumed by various sinks (memory, files, sharded outputs).
    """

    def __init__(self, source: PipelineSource) -> None:
        self._source: PipelineSource = source
        mgr = _lib.pipeline.StageManager(source)
        self._mgr = mgr
        # Let Python own the default policy - auto-inject built-ins only when referenced
        self._memory_sinks:     dict[str, list[_lib.pipeline.MemorySink]] = {}
        self._file_sinks:       dict[str, list[_lib.pipeline.NdjsonSink]] = {}
        self._sharded_sinks:    dict[str, list[_lib.pipeline.ShardedNdjsonSink]] = {}
        self._channel_formats:  dict[str, OutputFormat] = {}
        self._channel_decoders: dict[str, str] = {}

        # Parameter proxies
        self._system_params   = SystemParams(self._mgr)
        self._topology_params = TopologyParams(self._mgr)

        # hack
        self._contact_providers: set[_lib.ContactProvider] = set()

        # DatabaseHandleSource is only valid for is_model=True
        if isinstance(source, _lib.pipeline.sources.DatabaseHandleSource):
            try:
                self._system_params.is_model = True
            except Exception:
                pass
        self._reporting_level: ReportingLevel = ReportingLevel.BASIC
        self._last_report: dict[str, Any] | None = None

        try:
            self._mgr.set_auto_builtins(True)
        except Exception:
            pass
        try:
            self._mgr.set_flush_timeout(1.0) # To keep test runs fast
        except Exception:
            pass
        try:
            self._mgr.set_reporting_level(self._reporting_level)
        except AttributeError:
            # Older runtimes lack reporting level support. leave at default.
            # Check! Should not happen anymore.
            pass

    def set_reporting_level(self, level: ReportingLevel) -> None:
        """Control pipeline metrics collection."""
        if not isinstance(level, ReportingLevel):
            raise TypeError("level must be a ReportingLevel enum or string name")
        self._reporting_level = level
        self._mgr.set_reporting_level(level)

    def get_reporting_level(self) -> ReportingLevel:
        return self._reporting_level

    @staticmethod
    def _clone_report(report: dict[str, Any] | None) -> dict[str, Any] | None:
        if report is None:
            return None
        clone  = dict(report)
        stages = clone.get("stage_breakdown")
        if isinstance(stages, list):
            clone["stage_breakdown"] = [dict(stage) for stage in stages]
        return clone

    def get_run_report(self) -> dict[str, Any] | None:
        """
        Return the most recent run report, if available.

        When metrics are enabled the report includes pipeline totals as well as:
            * peak_inflight_items / average_queue_depth
            * permit_wait_* fields with min/avg/max acquisition latency
            * mux_* counters summarising channel writer activity
            * stage_breakdown (DEBUG): per-stage setup/compute timings
        """
        return self._clone_report(self._last_report)

    @overload
    def params(self, builtin_name: Literal["system"]) -> SystemParams: ...

    @overload
    def params(self, builtin_name: Literal["topology"]) -> TopologyParams: ...

    # Parameter access
    def params(self, builtin_name: str) -> SystemParams | TopologyParams:
        """Access parameters for built-in computations.

        Args:
            builtin_name: "system" or "topology"

        Returns:
            Parameter proxy object for the specified built-in

        Example:
            p.params("system").is_model = True
            p.params("topology").flags = TopologyComputation.Bonds
        """
        if builtin_name == "system":
            return self._system_params
        elif builtin_name == "topology":
            return self._topology_params
        else:
            raise ValueError(f"Unknown built-in: {builtin_name}. Valid options: 'system', 'topology'")

    def describe(self) -> str:
        """
        Return a description of the pipeline graph.
        Shows task dependencies and wheather nodes are built-in or user-defined.
        """
        graph = self._mgr.describe_graph()
        lines = ["Pipeline Graph:"]
        builtin_names = {"system", "topology"}
        for name, deps in graph:
            is_builtin = name in builtin_names
            node_type = "built-in" if is_builtin else "user"
            deps_str = " -> ".join(deps) if deps else "no deps"
            lines.append(f"  {name} ({node_type}): {deps_str}")
        return "\n".join(lines)

    def add_task(
        self,
        *,
        name: str,
        task: Callable[[PipelineContext], Any] | ContactTask | "CppTask",
        out: Iterable[FileOutput | ShardedOutput] | None = None,
        in_memory_policy: InMemoryPolicy = InMemoryPolicy.Keep,
        channel: str | None = None,
        store: bool | None = None,
        depends: list[str] | None = None,
        thread_safe: bool = True,
        writer_threads: int | None = None,
        requires_fields: Iterable[DataField] | None = None,
    ) -> None:
        ch = channel or name
        deps = list(depends) if depends is not None else []

        # Default behavior: AST to validate arity and infer builtins. If topology
        # is not found via AST, augment with a single probe run. Then union.
        if callable(task) and depends is None:
            ast_found: set[str] = set()
            try:
                ast_found = _ast_infer_builtins_for_callable(task)
            except ValueError as e:
                # Enforce arity==1 strictly with clear context
                raise ValueError(f"Invalid Python task signature for '{name}': {e}")

            runtime_found: list[str] = []
            if "topology" not in ast_found:
                runtime_found = _discover_builtins_for_callable(task)

            merged: set[str] = set(ast_found)
            merged.update(runtime_found)
            for d in merged:
                if d not in deps:
                    deps.append(d)

        # Guard: if topology is disabled, disallow tasks that require it
        if ("topology" in deps) and (not self._topology_params.enabled):
            raise ValueError(
                "Topology is disabled via Pipeline.params('topology').enabled = False. Cannot add tasks depending on 'topology'."
            )

        req_fields_tuple = self._normalize_requires_fields(requires_fields)
        if req_fields_tuple:
            self._validate_requires_fields()

        # Resolve ContactTask spec to a compute-backed contacts node
        if isinstance(task, ContactTask):
            if not self._topology_params.enabled:
                raise ValueError("Topology is disabled. Cannot add ContactsTask which requires 'topology'.")

            if isinstance(self._source, MdTrajectoriesSource) and self._contact_providers and task.provider not in self._contact_providers:
                existing = next(iter(self._contact_providers))
                raise ValueError(
                    f"Cannot add ContactTask with provider '{task.provider.name}' when MD source already has "
                    f"ContactTask with provider '{existing.name}'. MD trajectories share topology across frames, "
                    f"and different providers require incompatible atom typing. Use a single provider for a Pipeline run."
                )
            self._contact_providers.add(task.provider)

            match task.provider:
                case _lib.ContactProvider.Arpeggio:
                    self._topology_params.atom_typing_method = _lib.AtomTypingMethod.Arpeggio
                case _lib.ContactProvider.GetContacts:
                    self._topology_params.atom_typing_method = _lib.AtomTypingMethod.GetContacts
                case _:
                    self._topology_params.atom_typing_method = _lib.AtomTypingMethod.MolStar

            emission_fmt = task.fmt if isinstance(task.fmt, OutputFormat) else OutputFormat(task.fmt)
            sink_override = self._preferred_sink_format(out)
            if emission_fmt is OutputFormat.BINARY and sink_override is not None:
                emission_fmt = sink_override

            self._mgr.add_contacts(
                name,
                task.provider,
                task.interaction_type,
                ch,
                emission_fmt.value,
                bool(thread_safe),
            )

            self._channel_formats[ch]  = emission_fmt
            self._channel_decoders[ch] = "contacts"
            self._attach_sinks(ch, in_memory_policy, out, channel_format=emission_fmt, writer_threads=writer_threads)
            self._apply_requires_fields(name, req_fields_tuple)
            return

        # Python callable task
        if callable(task):
            do_store = True if store is None else bool(store)
            # Serialize Python callable by default unless user explicitly marks thread_safe
            serialize = not bool(thread_safe)
            self._mgr.add_python(name, deps, task, ch, serialize=serialize, store=do_store)

            # Auto-detect format from return type annotation and AST analysis
            task_format = infer_python_task_format(task)
            self._channel_formats[ch] = task_format
            self._attach_sinks(ch, in_memory_policy, out, channel_format=task_format, writer_threads=writer_threads)
            self._apply_requires_fields(name, req_fields_tuple)
            return

        # Raw C++ task
        if isinstance(task, _lib.pipeline.Task):
            self._mgr.add_task(name, deps, task, bool(thread_safe))
            self._channel_formats.setdefault(ch, OutputFormat.TEXT)
            self._attach_sinks(ch, in_memory_policy, out, channel_format=self._channel_formats[ch], writer_threads=writer_threads)
            self._apply_requires_fields(name, req_fields_tuple)
            return

        raise TypeError("task must be a callable, ContactTask spec, or pipeline.Task instance")

    def _normalize_requires_fields(self, fields: Iterable[DataField] | None) -> tuple[DataField, ...]:
        if not fields:
            return tuple()
        normalized: list[DataField] = []
        for entry in fields:
            if isinstance(entry, DataField):
                normalized.append(entry)
            else:
                raise TypeError("requires_fields entries must be DataField enums")
        return tuple(normalized)

    def _validate_requires_fields(self) -> None:
        is_model = bool(self._mgr.get_system_params().get("is_model", False))
        if not is_model:
            raise ValueError(
                "requires_fields is only supported for model inputs. "
                "Set pipeline.params('system').is_model = True before adding the task."
            )

    def _apply_requires_fields(self, task_name: str, fields: Sequence[DataField]) -> None:
        if not fields:
            return
        self._mgr.set_task_data_fields(task_name, list(fields))

    def _preferred_sink_format(self, out: Iterable[FileOutput | ShardedOutput] | None) -> OutputFormat | None:
        if not out:
            return None
        for o in out:
            if o.fmt == OutputFormat.BINARY:
                raise ValueError("Binary file output is not supported for file or sharded sinks.")
            if o.fmt == OutputFormat.JSON:
                return OutputFormat.JSON
            if o.fmt == OutputFormat.TEXT:
                return OutputFormat.TEXT
        return None

    def _make_backpressure_config(self, writer_threads: int | None):
        if writer_threads is None:
            return None
        if isinstance(writer_threads, bool) or not isinstance(writer_threads, Integral):
            raise TypeError("writer_threads must be a positive integer")
        value = int(writer_threads)
        if value <= 0:
            raise ValueError("writer_threads must be a positive integer")
        cfg = _lib.pipeline.get_default_backpressure_config()
        cfg.writer_threads = value
        cfg.validate()
        return cfg

    def _attach_sinks(
        self,
        channel: str,
        in_memory_policy: InMemoryPolicy,
        out: Iterable[FileOutput | ShardedOutput] | None,
        *,
        channel_format: OutputFormat,
        writer_threads: int | None = None,
    ) -> None:
        """Attach memory and file sinks for the given channel respecting output formats."""
        sink_cfg = self._make_backpressure_config(writer_threads)

        def connect(target_channel: str, sink_obj: Any) -> None:
            if sink_cfg is None:
                self._mgr.connect_sink(target_channel, sink_obj)
            else:
                self._mgr.connect_sink(target_channel, sink_obj, sink_cfg)

        if in_memory_policy == InMemoryPolicy.Keep:
            if channel not in self._memory_sinks:
                ms = _lib.pipeline.MemorySink()
                connect(channel, ms)
                self._memory_sinks[channel] = [ms]
            self._channel_formats.setdefault(channel, channel_format)

        if not out:
            return

        for o in out:
            if o.fmt == OutputFormat.BINARY:
                raise ValueError("File and sharded outputs do not support binary payloads.")
            if isinstance(o, FileOutput):
                sink = _lib.pipeline.NdjsonSink(str(o.path))
                connect(channel, sink)
                self._file_sinks.setdefault(channel, []).append(sink)
            elif isinstance(o, ShardedOutput):
                sink = _lib.pipeline.ShardedNdjsonSink(str(o.out_dir), int(o.shard_size))
                connect(channel, sink)
                self._sharded_sinks.setdefault(channel, []).append(sink)

    # Sinks
    def to_files(self, task_or_channel: str, *, path: str | Path, fmt: OutputFormat = OutputFormat.JSON, writer_threads: int | None = None) -> None:
        if not isinstance(fmt, OutputFormat):
            raise TypeError("fmt must be an OutputFormat enum value")
        if fmt == OutputFormat.BINARY:
            raise ValueError("Binary output is not supported for NdjsonSink.")
        sink_cfg = self._make_backpressure_config(writer_threads)
        sink = _lib.pipeline.NdjsonSink(str(path))
        if sink_cfg is None:
            self._mgr.connect_sink(task_or_channel, sink)
        else:
            self._mgr.connect_sink(task_or_channel, sink, sink_cfg)
        self._file_sinks.setdefault(task_or_channel, []).append(sink)
        self._channel_formats.setdefault(task_or_channel, fmt)
        return

    def to_memory(self, task_or_channel: str, *, writer_threads: int | None = None) -> None:
        sink_cfg = self._make_backpressure_config(writer_threads)
        if task_or_channel not in self._memory_sinks:
            ms = _lib.pipeline.MemorySink()
            if sink_cfg is None:
                self._mgr.connect_sink(task_or_channel, ms)
            else:
                self._mgr.connect_sink(task_or_channel, ms, sink_cfg)
            self._memory_sinks[task_or_channel] = [ms]

            #
            # If the channel already has a format (e.g., from a previous task), we keep it.
            # Otherwise default to JSON.
            #
            self._channel_formats.setdefault(task_or_channel, OutputFormat.JSON)
        return

    def to_sharded_files(
        self,
        task_or_channel: str,
        *,
        out_dir: str | Path,
        fmt: OutputFormat = OutputFormat.JSON,
        shard_size: int = 1000,
        writer_threads: int | None = None,
    ) -> None:
        if not isinstance(fmt, OutputFormat):
            raise TypeError("fmt must be an OutputFormat enum value")
        if fmt == OutputFormat.BINARY:
            raise ValueError("Binary output is not supported for sharded NDJSON sinks.")

        sink_cfg = self._make_backpressure_config(writer_threads)
        sink = _lib.pipeline.ShardedNdjsonSink(str(out_dir), int(shard_size))
        if sink_cfg is None:
            self._mgr.connect_sink(task_or_channel, sink)
        else:
            self._mgr.connect_sink(task_or_channel, sink, sink_cfg)
        self._sharded_sinks.setdefault(task_or_channel, []).append(sink)
        self._channel_formats.setdefault(task_or_channel, fmt)
        return

    def run(
        self,
        threads: int = 8,
    ) -> PipelineResult:
        self._last_report = None
        # Reset memory sinks so results do not accumulate across runs
        for sinks in self._memory_sinks.values():
            for s in sinks:
                try:
                    s.clear()
                except Exception:
                    pass

        def _collect_memory_states() -> dict[str, _ChannelState]:
            states: dict[str, _ChannelState] = {}
            for channel, sinks in self._memory_sinks.items():
                fmt = self._channel_formats.get(channel, OutputFormat.JSON)
                buffer: list[Any] = []
                for sink in sinks:
                    if fmt is OutputFormat.BINARY:
                        buffer.extend(sink.result_bytes())
                    else:
                        buffer.extend(sink.result())
                states[channel] = _ChannelState(
                    format=fmt,
                    decoder=self._channel_decoders.get(channel),
                    data=tuple(buffer),
                )
            return states

        _start = time.perf_counter()
        self._mgr.run(int(threads))
        _elapsed = time.perf_counter() - _start
        logging.info(f"pipeline.run finished in {_elapsed*1000:.1f} ms ({_elapsed:.3f} s)")
        report = self._mgr.last_run_report()
        self._last_report = self._clone_report(report)
        return PipelineResult(_collect_memory_states())

    #
    # If we want strong static typing of the run() output, we must accept a
    # TypedDict schema describing the channels and their element types,
    # and cast the runtime result to that schema for the type checker.
    # At runtime, this is identical to run(). This is the best that *I*
    # know how to do in Python's type system.    - Besian, September 2025
    #
    def run_typed(self, schema: Type[SchemaT], threads: int = 8) -> SchemaT:
        result = self.run(threads)
        return cast(SchemaT, dict(result.items()))

    def file_outputs(self) -> dict[str, list[str]]:
        res: dict[str, list[str]] = {}
        for ch, sinks in self._file_sinks.items():
            lst = res.setdefault(ch, [])
            for s in sinks:
                lst.append(s.file())
        for ch, sinks in self._sharded_sinks.items():
            lst = res.setdefault(ch, [])
            for s in sinks:
                lst.extend(s.files())
        return res


__all__ = ["Pipeline"]
