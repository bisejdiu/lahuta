from __future__ import annotations

import json as _json
from pathlib import Path
from typing import (
    Any,
    Callable,
    Iterable,
    Literal,
    Mapping,
    Protocol,
    Sequence,
    Type,
    TypeAlias,
    TypeVar,
    cast,
    overload,
)

from lahuta.lib import lahuta as _lib

from ._discovery import _ast_infer_builtins_for_callable, _discover_builtins_for_callable
from .params import SystemParams, TopologyParams
from .tasks import ContactTask
from .types import FileOutput, InMemoryPolicy, OutputFormat, PipelineContext, ShardedOutput

# fmt: off
StageManager = _lib.pipeline.StageManager
CppTask = _lib.pipeline.Task

JsonPrimitive = str | int | float | bool | None
JsonObject: TypeAlias = Mapping[str, "JsonValue"]
JsonArray:  TypeAlias = Sequence["JsonValue"]
JsonValue:  TypeAlias = JsonPrimitive | JsonObject | JsonArray

# Callable protocol for Python tasks that accept a PipelineContext and return any value
# that can be serialized into a single ndjson record. The return type is expressed
# by parameterizing PyTaskFn (e.g., PyTaskFn[dict[str, int]] or PyTaskFn[str]).
R_co = TypeVar("R_co", covariant=True)
SchemaT = TypeVar("SchemaT")


class PyTaskFn(Protocol[R_co]):
    def __call__(self, __ctx: PipelineContext) -> R_co: ...


class Pipeline:
    """Pipeline builder and executor.

    - Users add builtin tasks and Python tasks that accept PipelineContext
    - Each task emits to a channel. Sinks subscribe to channels
    - run() returns memory results keyed by channel. JSON channels auto-parse per item
    """

    def __init__(self, mgr: StageManager) -> None:
        self._mgr = mgr
        # Let Python own the default policy - auto-inject built-ins only when referenced
        try:
            self._mgr.set_auto_builtins(True)
        except Exception:
            pass

        self._json_memory_channels: set[str] = set()
        self._memory_sinks:  dict[str, list[_lib.pipeline.MemorySink]] = {}
        self._file_sinks:    dict[str, list[_lib.pipeline.NdjsonSink]] = {}
        self._sharded_sinks: dict[str, list[_lib.pipeline.ShardedNdjsonSink]] = {}

        # Parameter proxies
        self._system_params   = SystemParams(self._mgr)
        self._topology_params = TopologyParams(self._mgr)

    # Sources
    @staticmethod
    def from_directory(path: str | Path, ext: str = "", recursive: bool = True, batch: int = 200) -> "Pipeline":
        mgr = _lib.pipeline.StageManager.from_directory(str(path), ext, recursive, int(batch))
        return Pipeline(mgr)

    @staticmethod
    def from_files(files: str | Path | Sequence[str | Path]) -> "Pipeline":
        if isinstance(files, (str, Path)):
            paths = [str(files)]
        else:
            paths = [str(p) for p in files]
        mgr = _lib.pipeline.StageManager.from_files(paths)
        return Pipeline(mgr)

    @staticmethod
    def from_filelist(path: str | Path) -> "Pipeline":
        mgr = _lib.pipeline.StageManager.from_filelist(str(path))
        return Pipeline(mgr)

    @staticmethod
    def from_nmr_files(files: str | Path | Sequence[str | Path]) -> "Pipeline":
        """Create a pipeline that streams multi-model NMR structures."""

        if isinstance(files, (str, Path)):
            paths = [str(files)]
        else:
            paths = [str(p) for p in files]
        mgr = _lib.pipeline.StageManager.from_nmr_files(paths)
        return Pipeline(mgr)

    @staticmethod
    def from_database(path: str | Path, batch: int = 1024) -> "Pipeline":
        """Create a pipeline reading item keys from an LMDB database."""
        db = _lib.db.Database(str(path))
        return Pipeline.from_database_handle(db, int(batch))

    @staticmethod
    def from_database_handle(db: "_lib.db.Database", batch: int = 1024) -> "Pipeline":
        """Create a pipeline reading item keys from an LMDB database handle."""
        mgr = _lib.pipeline.StageManager.from_database_handle(db, int(batch))
        p = Pipeline(mgr)
        # SystemRead needs to run in model mode for LMDB-backed sessions
        try:
            p.params("system").is_model = True
        except Exception:
            pass
        return p

    @staticmethod
    def from_md_trajectories(trajectories: Iterable[Mapping[str, Any] | Sequence[Any]]) -> "Pipeline":
        """Create a pipeline that streams MD trajectories (structure + XTC)."""

        def _coerce_xtc_paths(value: Any) -> list[str]:
            if isinstance(value, (str, Path)):
                return [str(value)]
            return [str(v) for v in value]

        specs: list[dict[str, Any]] = []
        for entry in trajectories:
            if isinstance(entry, Mapping):
                structure_value = entry.get("structure", entry.get("gro"))
                if structure_value is None:
                    raise ValueError(
                        "trajectory mapping must include a 'structure' (or 'gro') key"
                    )
                structure = str(structure_value)
                xtc_value = entry.get("xtc", entry.get("xtcs"))
                if xtc_value is None:
                    raise ValueError("trajectory mapping must include 'xtc' or 'xtcs'")
                xtc_paths = _coerce_xtc_paths(xtc_value)
                session_id = entry.get("id", entry.get("session_id", structure))
                ident = str(session_id)
            else:
                seq = list(entry)
                if len(seq) < 2:
                    raise ValueError(
                        "trajectory tuple must provide structure and xtc paths"
                    )
                structure = str(seq[0])
                xtc_paths = _coerce_xtc_paths(seq[1])
                ident = str(seq[2]) if len(seq) > 2 else structure
            if not xtc_paths:
                raise ValueError("each trajectory requires at least one XTC path")
            specs.append({"structure": structure, "xtc": xtc_paths, "id": ident})

        mgr = _lib.pipeline.StageManager.from_md_trajectories(specs)
        return Pipeline(mgr)

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
            raise ValueError(
                f"Unknown built-in: {builtin_name}. Valid options: 'system', 'topology'"
            )

    def describe(self) -> str:
        """Return a description of the pipeline graph.

        Shows task dependencies and whether nodes are built-in or user-defined.
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
        depends: list[str] | None = None,
        channel: str | None = None,
        thread_safe: bool = True,
        store: bool | None = None,
        in_memory_policy: InMemoryPolicy = InMemoryPolicy.Keep,
        out: Iterable[FileOutput | ShardedOutput] | None = None,
    ) -> "Pipeline":
        ch = channel or name
        explicit_depends_provided = depends is not None
        deps = list(depends) if explicit_depends_provided else []

        # Default behavior: AST to validate arity and infer builtins. If topology
        # is not found via AST, augment with a single probe run. Then union.
        if callable(task) and not explicit_depends_provided:
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

            # Set atom typing method based on provider before adding contacts
            if task.provider == _lib.ContactProvider.Arpeggio:
                self._topology_params.atom_typing_method = _lib.AtomTypingMethod.Arpeggio
            else:  # MolStar
                self._topology_params.atom_typing_method = _lib.AtomTypingMethod.Molstar

            # JSON is the default out format
            self._mgr.add_contacts(
                name, deps, task.provider, task.interaction_type, ch, "json", bool(thread_safe)
            )
            self._attach_sinks(ch, in_memory_policy, out, json_channel=True)
            return self

        # Python callable task
        if callable(task):
            do_store = True if store is None else bool(store)
            # Serialize Python callable by default unless user explicitly marks thread_safe
            serialize = not bool(thread_safe)
            self._mgr.add_python(name, deps, task, ch, serialize=serialize, store=do_store)
            self._attach_sinks(ch, in_memory_policy, out, json_channel=True)
            return self

        # Raw C++ task
        if isinstance(task, _lib.pipeline.Task):
            self._mgr.add_task(name, deps, task, bool(thread_safe))
            self._attach_sinks(ch, in_memory_policy, out, json_channel=False)
            return self

        raise TypeError("task must be a callable, ContactTask spec, or pipeline.Task instance")

    def _attach_sinks(
        self,
        channel: str,
        in_memory_policy: InMemoryPolicy,
        out: Iterable[FileOutput | ShardedOutput] | None,
        json_channel: bool,
    ) -> None:
        """Simplified sink attachment logic."""
        # Memory sink only for text/JSON channels. Binary channels (e.g., LMDB payloads)
        # must not attach a MemorySink, as pybind11 decodes std::string to str (utf-8).
        if json_channel and in_memory_policy == InMemoryPolicy.Keep:
            if channel not in self._memory_sinks:
                ms = _lib.pipeline.MemorySink()
                self._mgr.connect_sink(channel, ms)
                self._memory_sinks[channel] = [ms]
                self._json_memory_channels.add(channel)

        if not out:
            return

        # File/sharded sinks
        for o in out:
            if isinstance(o, FileOutput):
                sink = _lib.pipeline.NdjsonSink(str(o.path))
                self._mgr.connect_sink(channel, sink)
                self._file_sinks.setdefault(channel, []).append(sink)
            elif isinstance(o, ShardedOutput):
                sink = _lib.pipeline.ShardedNdjsonSink(str(o.out_dir), int(o.shard_size))
                self._mgr.connect_sink(channel, sink)
                self._sharded_sinks.setdefault(channel, []).append(sink)

    # Sinks
    def to_files(self, task_or_channel: str, *, path: str | Path, fmt: OutputFormat = OutputFormat.JSON) -> "Pipeline":
        if not isinstance(fmt, OutputFormat):
            raise TypeError("fmt must be an OutputFormat enum value")
        sink = _lib.pipeline.NdjsonSink(str(path))
        self._mgr.connect_sink(task_or_channel, sink)
        self._file_sinks.setdefault(task_or_channel, []).append(sink)
        return self

    def to_memory(self, task_or_channel: str) -> "Pipeline":
        if task_or_channel not in self._memory_sinks:
            ms = _lib.pipeline.MemorySink()
            self._mgr.connect_sink(task_or_channel, ms)
            self._memory_sinks[task_or_channel] = [ms]
        return self

    def to_sharded_files(
        self,
        task_or_channel: str,
        *,
        out_dir: str | Path,
        fmt: OutputFormat = OutputFormat.JSON,
        shard_size: int = 1000,
    ) -> "Pipeline":
        if not isinstance(fmt, OutputFormat):
            raise TypeError("fmt must be an OutputFormat enum value")
        sink = _lib.pipeline.ShardedNdjsonSink(str(out_dir), int(shard_size))
        self._mgr.connect_sink(task_or_channel, sink)
        self._sharded_sinks.setdefault(task_or_channel, []).append(sink)
        return self

    def run(self, threads: int = 8) -> dict[str, list[Any]]:
        #
        # We use Any to keep static analyzers permissive across heterogeneous channels.
        # At runtime, json channels are parsed to native Python values. Text channels
        # yield raw strings. Per-channel precise typing can be modeled by users py
        # parameterizing PyTaskFn on the callables they define.   - Besian, September 2025
        #

        # Reset memory sinks so results do not accumulate across runs
        for sinks in self._memory_sinks.values():
            for s in sinks:
                try:
                    s.clear()
                except Exception:
                    # fall back if clear() is not implemented
                    pass

        self._mgr.run(int(threads))
        out: dict[str, list[Any]] = {}
        for ch, sinks in self._memory_sinks.items():
            # We use Any to keep static analyzers permissive across heterogeneous channels.
            # At runtime, json channels are parsed to native Python values. Text channels
            # yield raw strings.
            acc: list[Any] = []
            for s in sinks:
                vals = s.result()
                if ch in self._json_memory_channels:
                    for v in vals:
                        try:
                            acc.append(_json.loads(v))
                        except Exception:
                            acc.append(v)
                else:
                    acc.extend(vals)
            out[ch] = acc
        return out

    #
    # If we want strong static typing of the run() output, we must accept a
    # TypedDict schema describing the channels and their element types,
    # and cast the runtime result to that schema for the type checker.
    # At runtime, this is identical to run(). This is the best that *I*
    # know how to do in Python's type system.    - Besian, September 2025
    #
    def run_typed(self, schema: Type[SchemaT], threads: int = 8) -> SchemaT:
        return cast(SchemaT, self.run(threads))

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
