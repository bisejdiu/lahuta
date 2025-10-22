from __future__ import annotations

import os
import time
from numbers import Integral
from pathlib import Path
from typing import (
    Any,
    Callable,
    Iterable,
    Literal,
    Protocol,
    Type,
    TypeAlias,
    TypeVar,
    cast,
    overload,
)

from lahuta import logging
from lahuta.lib import lahuta as _lib
from lahuta.sources import PipelineSource

from ._discovery import _ast_infer_builtins_for_callable, _discover_builtins_for_callable
from ._inference import infer_python_task_format
from .mp_backend import PyTaskSpec as _MPPyTaskSpec
from .mp_backend import resolve_callable_metadata
from .params import SystemParams, TopologyParams
from .result import PipelineResult, _ChannelState
from .tasks import ContactTask
from .types import FileOutput, InMemoryPolicy, OutputFormat, PipelineContext, ShardedOutput

# fmt: off
StageManager = _lib.pipeline.StageManager
ReportingLevel = _lib.pipeline.ReportingLevel
CppTask: TypeAlias = _lib.pipeline.Task

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

        # DatabaseHandleSource is only valid for is_model=True
        if isinstance(source, _lib.pipeline.sources.DatabaseHandleSource):
            try:
                self._system_params.is_model = True
            except Exception:
                pass
        self._py_tasks: list[_MPPyTaskSpec] = [] # callable task registry for multiprocessing backend
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
            # Older runtimes lack reporting level support; leave at default.
            pass

    def set_reporting_level(self, level: ReportingLevel | str) -> None:
        """Control pipeline metrics collection."""
        if isinstance(level, str):
            try:
                level = ReportingLevel[level.upper()]
            except KeyError as exc:
                raise ValueError(f"Unknown reporting level '{level}'. Valid levels: {[e.name for e in ReportingLevel]}") from exc
        if not isinstance(level, ReportingLevel):
            raise TypeError("level must be a ReportingLevel enum or string name")
        self._reporting_level = level
        try:
            self._mgr.set_reporting_level(level)
        except AttributeError:
            # Older extensions do not expose reporting levels
            raise RuntimeError("Current Lahuta extension does not support reporting levels; rebuild Lahuta.") from None

    def get_reporting_level(self) -> ReportingLevel:
        return self._reporting_level

    def get_run_report(self) -> dict[str, Any] | None:
        """Return the most recent run report, if available."""
        return None if self._last_report is None else dict(self._last_report)

    def _process_pool_guard(self, processes: int):
        class _Guard:
            def __init__(self, mgr, count):
                self._mgr   = mgr
                self._count = count

            def __enter__(self):
                self._mgr.configure_python_process_pool(self._count)
                return self

            def __exit__(self, exc_type, exc, tb):
                try:
                    self._mgr.shutdown_python_process_pool()
                except Exception:
                    pass

        return _Guard(self._mgr, processes)

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

        # Resolve ContactTask spec to a compute-backed contacts node
        if isinstance(task, ContactTask):
            if not self._topology_params.enabled:
                raise ValueError("Topology is disabled. Cannot add ContactsTask which requires 'topology'.")

            match task.provider:
                case _lib.ContactProvider.Arpeggio:
                    self._topology_params.atom_typing_method = _lib.AtomTypingMethod.Arpeggio
                case _lib.ContactProvider.GetContacts:
                    self._topology_params.atom_typing_method = _lib.AtomTypingMethod.GetContacts
                case _:
                    self._topology_params.atom_typing_method = _lib.AtomTypingMethod.Molstar

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

            # Record for multiprocessing backend
            module, qualname = resolve_callable_metadata(task)
            callable_blob: bytes | None = None
            if module is None or qualname is None:
                # Notebook/REPL task - serialize via cloudpickle
                try:
                    import cloudpickle
                    callable_blob = cloudpickle.dumps(task)
                except ImportError:
                    logging.warn(
                        f"Process backend: cloudpickle not available, task '{name}' will only run in threaded mode. "
                        "Install cloudpickle to enable multiprocessing support for notebook tasks."
                    )
                    callable_blob = None
                except Exception as exc:
                    logging.warn(
                        f"Process backend: unable to serialize task '{name}' via cloudpickle ({exc}). "
                        "Task will only run in threaded mode."
                    )
                    callable_blob = None

            self._py_tasks.append(_MPPyTaskSpec(
                name=name,
                fn=cast(PyTaskFn[Any], task),
                module=module,
                qualname=qualname,
                callable_blob=callable_blob,
                channel=ch,
                store=do_store,
                serialize=serialize,
                depends=tuple(deps),
                format=task_format,
            ))
            return

        # Raw C++ task
        if isinstance(task, _lib.pipeline.Task):
            self._mgr.add_task(name, deps, task, bool(thread_safe))
            self._channel_formats.setdefault(ch, OutputFormat.TEXT)
            self._attach_sinks(ch, in_memory_policy, out, channel_format=self._channel_formats[ch], writer_threads=writer_threads)
            return

        raise TypeError("task must be a callable, ContactTask spec, or pipeline.Task instance")

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
        *,
        backend: Literal["threads", "processes"] = "threads",
        processes: int | None = None,
        process_timeout: float | None = 300.0,
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

        if backend == "threads":
            _start = time.perf_counter()
            self._mgr.run(int(threads))
            _elapsed = time.perf_counter() - _start
            logging.info(f"pipeline.run finished in {_elapsed*1000:.1f} ms ({_elapsed:.3f} s)")
            report = self._mgr.last_run_report()
            self._last_report = dict(report) if report is not None else None
            return PipelineResult(_collect_memory_states())

        if backend == "processes":
            if not self._py_tasks:
                logging.warn("No Python tasks registered. Falling back to threaded backend.")
                return self.run(threads=threads, backend="threads")

            streaming_types: tuple[type[Any], ...] = ()
            try:
                # Fast failure until process backend supports streaming/MD descriptors across processes.
                streaming_types = _lib.pipeline.sources.MdTrajectoriesSource, _lib.pipeline.sources.NmrSource
            except AttributeError:
                streaming_types = ()

            if streaming_types and isinstance(self._source, streaming_types):
                raise RuntimeError("backend='processes' is not supported for streaming MD/NMR sources.")

            # Tasks are convertible if they have module/qualname OR serialized blob
            convertible = [
                spec for spec in self._py_tasks
                if (spec.module and spec.qualname) or (spec.callable_blob is not None)
            ]
            if not convertible:
                logging.warn("No process-capable Python tasks registered. Falling back to threaded backend.")
                return self.run(threads=threads, backend="threads")

            proc_count      = processes if processes is not None else os.cpu_count() or 1
            proc_count      = max(1, int(proc_count))
            worker_threads  = max(1, int(threads))
            timeout_seconds = float(process_timeout) if process_timeout is not None else 300.0

            converted_names: set[str] = set()
            try:
                with self._process_pool_guard(proc_count):
                    self._mgr.set_python_process_concurrency(proc_count)
                    self._mgr.set_python_process_timeout(timeout_seconds)
                    for spec in convertible:
                        # Pass serialized callable as bytes or None
                        serialized = spec.callable_blob if spec.callable_blob is not None else None
                        self._mgr.add_python_process(
                            spec.name,
                            list(spec.depends),
                            spec.module or "",
                            spec.qualname or "",
                            serialized,
                            spec.channel,
                            store=spec.store,
                        )
                        converted_names.add(spec.name)

                    _start = time.perf_counter()
                    self._mgr.run(worker_threads)
                    _elapsed = time.perf_counter() - _start
                    logging.info(f"pipeline.run (processes) finished in {_elapsed*1000:.1f} ms ({_elapsed:.3f} s)")
            finally:
                if converted_names:
                    for spec in self._py_tasks:
                        if spec.name in converted_names:
                            self._mgr.add_python(
                                spec.name,
                                list(spec.depends),
                                spec.fn,
                                spec.channel,
                                serialize=spec.serialize,
                                store=spec.store,
                            )

            report = self._mgr.last_run_report()
            self._last_report = dict(report) if report is not None else None
            return PipelineResult(_collect_memory_states())

        raise ValueError("backend must be either 'threads' or 'processes'")

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
