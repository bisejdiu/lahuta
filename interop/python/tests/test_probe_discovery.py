# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     d = collections.deque()
#     d.appendleft("@gmail.com")
#     d.appendleft("sejdiu")
#     d.appendleft("besian")
#     print("".join(d))
#
from __future__ import annotations

from lahuta.pipeline._discovery import _discover_builtins_for_callable
from lahuta.pipeline._probe import (
    _ProbePipelineContext,
    _ProbeValue,
)


def test_probe_value_truthiness_and_operations() -> None:
    v = _ProbeValue("X")
    # False-y by default
    assert bool(v) is False
    assert str(v) == ""
    assert int(v) == 0
    assert float(v) == 0.0
    # Iteration is empty
    assert list(v) == []
    # Arithmetic returns the same probe value
    assert (v + 1) is v
    assert (1 + v) is v
    assert (v * 2) is v
    assert (v / 2) is v


def test_probe_context_records_accessors_direct_calls() -> None:
    ctx = _ProbePipelineContext()
    s = ctx.get_system()
    t = ctx.get_topology()
    assert bool(s) is False
    assert bool(t) is False
    assert ctx._used == {"system", "topology"}


def test_probe_context_records_accessors_via_getattr() -> None:
    ctx = _ProbePipelineContext()
    getattr(ctx, "get_system")()
    getattr(ctx, "get_topology")()
    assert ctx._used == {"system", "topology"}


def test_probe_value_method_calls_and_chaining() -> None:
    ctx = _ProbePipelineContext()
    sys = ctx.get_system()
    # Method calls return the same probe value, do not raise
    n = sys.get_number_of_atoms()
    res = sys.run_computation()
    chain = sys.foo.bar().baz(1, key="x")
    # All are probe values and falsey
    assert isinstance(n, _ProbeValue)
    assert isinstance(res, _ProbeValue)
    assert isinstance(chain, _ProbeValue)
    assert not n and not res and not chain


def test_discovery_returns_system_only_when_only_system_used() -> None:
    def needs_system(ctx: _ProbePipelineContext):
        _ = ctx.get_system()
        return None

    deps = _discover_builtins_for_callable(needs_system)
    assert deps == ["system"]


def test_discovery_returns_topology_when_topology_used() -> None:
    def needs_topology(ctx: _ProbePipelineContext):
        _ = ctx.get_topology()
        return None

    deps = _discover_builtins_for_callable(needs_topology)
    assert deps == ["topology"]


def test_discovery_ignores_functions_with_no_accessors() -> None:
    def no_access(_ctx: _ProbePipelineContext):
        return "ok"

    deps = _discover_builtins_for_callable(no_access)
    assert deps == []


def test_discovery_captures_before_exception() -> None:
    def raises_after_access(ctx: _ProbePipelineContext):
        _ = ctx.get_system()
        raise RuntimeError("boom")

    deps = _discover_builtins_for_callable(raises_after_access)
    assert deps == ["system"]


def test_discovery_misses_untaken_branch_expected() -> None:
    # Probe-only discovery will not see access hidden behind an untaken branch
    def hidden(x: _ProbePipelineContext):
        if False:
            _ = x.get_topology()
        return None

    deps = _discover_builtins_for_callable(hidden)
    assert deps == []  # Expected miss. We cannot infer topology dependency as the probe does not execute that branch.


def test_discovery_dynamic_getattr_name_expr() -> None:
    # Name built dynamically. Probe still executes and records access
    def dyn(ctx: _ProbePipelineContext):
        name = "get_" + "system"
        getattr(ctx, name)()
        return None

    deps = _discover_builtins_for_callable(dyn)
    assert deps == ["system"]
