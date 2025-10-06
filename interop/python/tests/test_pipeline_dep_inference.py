from __future__ import annotations

from pathlib import Path
from typing import cast

import pytest

from lahuta.pipeline import Pipeline, PipelineContext, PyTaskFn
from lahuta.sources import FilesSource

DATA_FILE = Path(__file__).resolve().parents[3] / "core" / "data" / "1kx2_small.cif"


def _extract_task_deps(p: Pipeline, task_name: str) -> list[str]:
    # Describe graph returns [(name, [deps...]), ...]
    graph = p._mgr.describe_graph()
    for name, deps in graph:
        if name == task_name:
            return list(deps)
    raise AssertionError(f"Task {task_name!r} not found in graph: {graph}")


def _pipeline_for_test() -> Pipeline:
    return Pipeline(FilesSource([str(DATA_FILE)]))


def test_infer_system_direct_call() -> None:
    def f(ctx: PipelineContext):
        _ = ctx.get_system()
        return {"ok": True}

    p = _pipeline_for_test()
    p.add_task(name="f", task=f)
    assert _extract_task_deps(p, "f") == ["system"]


def test_infer_topology_direct_call() -> None:
    def g(ctx: PipelineContext):
        _ = ctx.get_topology()
        return {"ok": True}

    p = _pipeline_for_test()
    p.add_task(name="g", task=g)
    assert _extract_task_deps(p, "g") == ["topology"]


def test_infer_via_getattr_literal() -> None:
    def h(ctx: PipelineContext):
        getattr(ctx, "get_topology")()
        return None

    p = _pipeline_for_test()
    p.add_task(name="h", task=h)
    assert _extract_task_deps(p, "h") == ["topology"]


def test_infer_with_alias_assignment() -> None:
    def k(ctx: PipelineContext):
        c = ctx
        c2 = c
        _ = c2.get_system()
        return None

    p = _pipeline_for_test()
    p.add_task(name="k", task=k)
    assert _extract_task_deps(p, "k") == ["system"]


def test_infer_from_callable_class_dunder_call() -> None:
    class C:
        def __call__(self, ctx: PipelineContext):
            # Should infer topology
            return bool(ctx.get_topology())

    p = _pipeline_for_test()
    p.add_task(name="c", task=C())
    assert _extract_task_deps(p, "c") == ["topology"]


def test_infer_from_type_dynamic_class() -> None:
    # type-based dynamic class with __call__ (self bound -> one visible param)
    def __call__(self, ctx: PipelineContext):
        return bool(ctx.get_system())

    T = type("T", (), {"__call__": __call__})

    p = _pipeline_for_test()
    p.add_task(name="t", task=cast(PyTaskFn, T()))
    assert _extract_task_deps(p, "t") == ["system"]


def test_infer_runtime_only_dynamic_name_concat() -> None:
    def r(ctx: PipelineContext):
        name = "get_" + "system"
        getattr(ctx, name)()
        return None

    p = _pipeline_for_test()
    p.add_task(name="r", task=r)
    # AST will not see constant name; runtime probe should catch 'system'
    assert _extract_task_deps(p, "r") == ["system"]


def test_infer_ast_covers_untaken_branch() -> None:
    def s(ctx: PipelineContext):
        if False:
            _ = ctx.get_topology()
        return None

    p = _pipeline_for_test()
    p.add_task(name="s", task=s)
    # AST should still find the topology accessor
    assert _extract_task_deps(p, "s") == ["topology"]


def test_explicit_depends_overrides_inference() -> None:
    def z(ctx: PipelineContext):
        _ = ctx.get_topology()
        return None

    p = _pipeline_for_test()
    # Explicitly restrict to system only
    p.add_task(name="z", task=z, depends=["system"])
    assert _extract_task_deps(p, "z") == ["system"]


def test_arity_violation_raises() -> None:
    def bad(a, b):
        return None

    p = _pipeline_for_test()
    with pytest.raises(ValueError):
        p.add_task(name="bad", task=bad)  # type: ignore # ignore because we intentionally violate the signature
