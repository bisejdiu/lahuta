from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

from .types import OutputFormat


@dataclass(frozen=True)
class PyTaskSpec:
    """Metadata describing a Python task registered with the pipeline.

    The callable itself is retained so the threaded backend can continue to
    execute tasks in-process, while module/qualname are propagated to the
    C++ multiprocessing implementation to resolve the callable inside worker
    processes. If a callable cannot be referenced by module/qualname (common
    for notebook or REPL definitions), `callable_blob` holds a cloudpickle
    payload that the worker can deserialize instead.
    """

    name: str
    fn: Callable[[Any], Any]
    module: str | None
    qualname: str | None
    callable_blob: bytes | None
    channel: str
    store: bool
    serialize: bool
    depends: tuple[str, ...]
    format: OutputFormat


def resolve_callable_metadata(fn: Callable[[Any], Any]) -> tuple[str | None, str | None]:
    """Return (module, qualname) for a Python callable suitable for import.

    Returns (None, None) when the callable cannot be resolved to an importable
    symbol (e.g., nested functions, callables synthesized via `type(...)`).
    """
    module = getattr(fn, "__module__", None)
    qualname = getattr(fn, "__qualname__", None)
    if not module or not qualname:
        return None, None
    # Reject callables defined in __main__, notebooks, or REPL contexts
    if module == "__main__":
        return None, None
    if module.startswith("ipykernel") or module.startswith("IPython."):
        return None, None
    if module.startswith("<"):
        return None, None
    if "<locals>" in qualname:
        return None, None
    return module, qualname


__all__ = ["PyTaskSpec", "resolve_callable_metadata"]
