from __future__ import annotations

from typing import Any as _Any


class _ProbeValue:
    def __init__(self, label: str = "PROBE") -> None:
        self._label = label

    def __repr__(self) -> str:
        return f"<{self._label}>"

    def __str__(self) -> str:
        return ""

    def __bool__(self) -> bool:
        return False

    def __int__(self) -> int:
        return 0

    def __float__(self) -> float:
        return 0.0

    def __len__(self) -> int:
        return 0

    def __iter__(self):
        return iter(())

    def __call__(self, *a: _Any, **k: _Any):
        return self

    def __getattr__(self, _name: str):
        return self

    def _op(self, *_, **__):
        return self

    __add__ = __radd__ = __sub__ = __rsub__ = _op
    __mul__ = __rmul__ = __truediv__ = __rtruediv__ = _op
    __floordiv__ = __rfloordiv__ = __mod__ = __rmod__ = _op
    __pow__ = __rpow__ = _op


class _ProbePipelineContext:
    """Lightweight context used only for one-shot dependency discovery.

    Records which accessors were called. Returns permissive probe values to
    let user callables run without touching C++ or external state.
    """

    def __init__(self) -> None:
        self._used: set[str] = set()
        self.path = "__probe__"

    def get_system(self):
        self._used.add("system")
        return _ProbeValue("SYSTEM")

    def get_topology(self):
        self._used.add("topology")
        return _ProbeValue("TOPOLOGY")

    def get_text(self, _key: str):
        return _ProbeValue("TEXT")

    def get_json(self, _key: str):
        return _ProbeValue("JSON")

    def get(self, _key: str):
        return _ProbeValue("VAL")

    def set_text(self, _key: str, _val: str) -> None:
        return None

    def set_json(self, _key: str, _val: _Any) -> None:
        return None


__all__ = ["_ProbeValue", "_ProbePipelineContext"]
