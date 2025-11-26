from __future__ import annotations

import importlib
from typing import Any, Callable

import orjson

from lahuta.lib import lahuta as _lib

# fmt: off
_CALLABLE_CACHE: dict[tuple[str, str], Callable[..., Any]] = {}


def _resolve_callable(module: str, qualname: str) -> Callable[..., Any]:
    key = (module, qualname)
    fn = _CALLABLE_CACHE.get(key)
    if fn is not None:
        return fn

    target: Any = importlib.import_module(module)
    for part in qualname.split("."):
        target = getattr(target, part)
    result: Callable[..., Any] = target
    _CALLABLE_CACHE[key] = result
    return result


class _WorkerPipelineContext:
    """Lightweight substitute for PipelineContext usable inside worker processes."""

    def __init__(
        self,
        path: str,
        system: _lib.LahutaSystem,
        topology: _lib.Topology,
        initial_store: dict[str, str] | None,
        initial_bytes: dict[str, bytes] | None,
    ) -> None:
        self._path = path
        self._system   = system
        self._topology = topology
        self._store:  dict[str, str] = dict(initial_store or {})
        self._bytes:  dict[str, bytes] = dict(initial_bytes or {})
        self._writes: dict[str, str] = {}
        self._writes_bytes: dict[str, bytes] = {}

    @property
    def path(self) -> str:
        return self._path

    def conformer_id(self) -> int:
        return 0

    def session_id(self) -> None:
        return None

    def timestamp_ps(self) -> None:
        return None

    def frame_metadata(self) -> None:
        return None

    def get_system(self) -> _lib.LahutaSystem:
        return self._system

    def get_topology(self) -> _lib.Topology:
        return self._topology

    def set_text(self, key: str, value: str) -> None:
        self._store[key]  = value
        self._writes[key] = value

    def set_json(self, key: str, value: Any) -> None:
        self.set_text(key, orjson.dumps(value, option=orjson.OPT_SERIALIZE_NUMPY).decode("utf-8"))

    def get_text(self, key: str) -> str | None:
        return self._store.get(key)

    def get_json(self, key: str) -> Any | None:
        raw = self._store.get(key)
        if raw is None:
            return None
        try:
            return orjson.loads(raw)
        except Exception:
            return None

    def get(self, key: str) -> Any | None:
        val = self.get_json(key)
        if val is not None:
            return val
        return self.get_text(key)

    def set_bytes(self, key: str, value: bytes) -> None:
        self._bytes[key] = value
        self._writes_bytes[key] = value

    def get_bytes(self, key: str) -> bytes | None:
        return self._bytes.get(key)

    def flush_writes(self) -> dict[str, str]:
        writes = dict(self._writes)
        self._writes.clear()
        return writes

    def flush_writes_bytes(self) -> dict[str, bytes]:
        writes = dict(self._writes_bytes)
        self._writes_bytes.clear()
        return writes


_BASE_TOPOLOGY_FLAGS = (
    _lib.TopologyComputers.Neighbors,
    _lib.TopologyComputers.Bonds,
    _lib.TopologyComputers.NonStandardBonds,
    _lib.TopologyComputers.Residues,
    _lib.TopologyComputers.Rings,
    _lib.TopologyComputers.AtomTyping,
)


def _as_optional_atom_typing(method: Any) -> _lib.AtomTypingMethod | None:
    if method is None:
        return None
    try:
        value = int(method)
    except (TypeError, ValueError):
        return None
    if value < 0:
        return None
    try:
        return _lib.AtomTypingMethod(value)
    except Exception:
        return None


def _as_optional_flags(flags: Any) -> int | None:
    if flags is None:
        return None
    try:
        value = int(flags)
    except (TypeError, ValueError):
        return None
    if value < 0:
        return None
    return value

#
# Workers always rebuild topology because LMDB-backed DataField views
# cannot be shared across processes. Once we support copying the requested
# slices before dispatch, we can skip this rebuild. - Besian, November 2025
#

def _create_system(item_path: str, system_params: dict[str, Any] | None) -> _lib.LahutaSystem:
    is_model = bool(system_params.get("is_model")) if system_params else False
    if is_model:
        return _lib.LahutaSystem.from_model_file(item_path)
    return _lib.LahutaSystem(item_path)


def _apply_topology_params(
    system: _lib.LahutaSystem,
    topology_params: dict[str, Any] | None,
    is_model_mode: bool,
) -> bool:
    opts = _lib.TopologyBuildingOptions()
    flags_value = _as_optional_flags(topology_params.get("flags")) if topology_params else None
    atom_typing = _as_optional_atom_typing(topology_params.get("atom_typing_method") if topology_params else None)

    if atom_typing is not None:
        opts.atom_typing_method = atom_typing
        system.set_atom_typing_method(atom_typing)

    if is_model_mode and hasattr(_lib, "TopologyBuildMode"):
        try:
            opts.mode = _lib.TopologyBuildMode.Model  # type: ignore[attr-defined]
        except Exception:
            pass

    if flags_value is not None:
        opts.compute_nonstandard_bonds = bool(flags_value & int(_lib.TopologyComputers.NonStandardBonds))

    if not system.build_topology(opts):
        return False

    if flags_value is not None:
        try:
            system.enable_only(_lib.TopologyComputers(flags_value))
        except Exception:
            # Fallback: enable each base flag explicitly
            mask = 0
            for flag in _BASE_TOPOLOGY_FLAGS:
                if flags_value & int(flag):
                    mask |= int(flag)
            system.enable_only(_lib.TopologyComputers(mask))
    return True


def execute_callable(
    module: str,
    qualname: str,
    serialized_callable: bytes | None,
    item_path: str,
    initial_store: dict[str, str] | None,
    initial_bytes: dict[str, bytes] | None,
    system_params: dict[str, Any] | None = None,
    topology_params: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Worker entry point invoked from the C++ multiprocessing bridge."""
    ctx: _WorkerPipelineContext | None = None
    try:
        # Resolve callable: prefer serialized blob for notebook/REPL tasks
        if serialized_callable:
            import cloudpickle
            fn = cloudpickle.loads(serialized_callable)
        else:
            fn = _resolve_callable(module, qualname)
        system = _create_system(item_path, system_params)
        is_model_flag = bool(system_params.get("is_model")) if system_params else False
        if not _apply_topology_params(system, topology_params, is_model_flag):
            return {
                "ok": False,
                "payload": _serialize_error_payload(item_path, "Failed to build topology"),
                "store_text": {},
            }

        topology = system.get_topology()
        if topology is None:
            return {
                "ok": False,
                "payload": _serialize_error_payload(item_path, "Topology is None after build"),
                "store_text": {},
            }

        ctx = _WorkerPipelineContext(item_path, system, topology, initial_store, initial_bytes)
        result = fn(ctx)

        # Build typed response
        payload_text:  str | None
        payload_bytes: bytes | None
        if result is None:
            payload_text  = None
            payload_bytes = None
        elif isinstance(result, str):
            payload_text  = result
            payload_bytes = None
        elif isinstance(result, (bytes, bytearray, memoryview)):
            payload_text  = None
            payload_bytes = bytes(result)
        else:
            payload_text = orjson.dumps(result, option=orjson.OPT_SERIALIZE_NUMPY).decode("utf-8")
            payload_bytes = None

        response: dict[str, Any] = {
            "version": 1,
            "ok": True,
            "store_text": ctx.flush_writes(),
            "store_bytes": ctx.flush_writes_bytes(),
        }
        if payload_bytes is not None:
            response["payload_bytes"] = payload_bytes
        else:
            response["payload"] = payload_text
        return response
    except Exception as exc:
        message = f"{type(exc).__name__}: {str(exc).splitlines()[0]}"
        payload = _serialize_error_payload(item_path, message)
        writes = ctx.flush_writes() if ctx is not None else {}
        return {
            "version": 1,
            "ok": False,
            "payload": payload,
            "store_text": writes,
            "store_bytes": {},
        }


def _serialize_error_payload(item_path: str, message: str) -> str:
    payload = {
        "error": {"source": "python", "message": message},
        "path": item_path,
    }
    return orjson.dumps(payload, option=orjson.OPT_SERIALIZE_NUMPY).decode("utf-8")


__all__ = ["execute_callable"]
