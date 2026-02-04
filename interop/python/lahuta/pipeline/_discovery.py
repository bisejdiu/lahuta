# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     class Email:
#         def __init__(self, v): self.v = v
#         def __or__(self, f): return Email(f(self.v))
#         def __str__(self): return self.v
#     print(str(
#         Email("")
#         | (lambda s: s + "besian")
#         | (lambda s: s + "sejdiu")
#         | (lambda s: s + "@gmail.com")
#     ))
#
from __future__ import annotations

import ast
import inspect
from typing import Any

from ._probe import _ProbePipelineContext


# fmt: off
def _discover_builtins_for_callable(func: Any) -> list[str]:
    """Run func once with a probe context and return inferred built-in deps.

    Rules:
    - If topology is used, return ["topology"]
    - Else if system is used, return ["system"].
    - Else, return [].
    """
    probe = _ProbePipelineContext()
    try:
        func(probe)
    except Exception:
        pass
    if "topology" in probe._used: return ["topology"]
    if "system"   in probe._used: return ["system"]
    return []


def _ast_infer_builtins_for_callable(func: Any) -> set[str]:
    """Infer built-in dependencies

    Rules/limits:
    - Validates that the callable has exactly one non-vararg parameter.
    - Tracks the parameter name and simple aliases: a = ctx; b = a
    - Detects calls of the form ctx.get_system()/get_topology()
      and getattr(ctx, "get_system"/"get_topology")()
    - Returns a set containing any of {"system", "topology"}
    """
    sig = inspect.signature(func)
    params = [p for p in sig.parameters.values() if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)]
    if len(params) != 1 or any(p.kind in (p.VAR_POSITIONAL, p.VAR_KEYWORD) for p in sig.parameters.values()):
        raise ValueError(
            "Python task must accept exactly one positional parameter (PipelineContext). "
            f"Got signature: {func.__name__}{sig}"
        )
    arg_name = params[0].name

    try:
        src = inspect.getsource(func)
    except Exception:
        return set()

    try:
        tree = ast.parse("if 1:\n" + src) # textwrap.dedent(src) alsow works fine for ast.parse
    except Exception:
        return set()

    # Track aliases of the context argument via simple assignments
    aliases: set[str] = {arg_name}

    class AliasCollector(ast.NodeVisitor):
        def visit_Assign(self, node: ast.Assign) -> None:
            try:
                if isinstance(node.value, ast.Name) and node.value.id in aliases:
                    for tgt in node.targets:
                        if isinstance(tgt, ast.Name):
                            aliases.add(tgt.id)
            finally:
                self.generic_visit(node)

    class DepCollector(ast.NodeVisitor):
        def __init__(self) -> None:
            self.found: set[str] = set()

        def visit_Call(self, node: ast.Call) -> None:
            # Attribute form: <alias>.get_system()/get_topology()
            f = node.func
            if isinstance(f, ast.Attribute) and isinstance(f.value, ast.Name) and f.value.id in aliases:
                if f.attr == "get_system":
                    self.found.add("system")
                elif f.attr == "get_topology":
                    self.found.add("topology")

            # getattr form: getattr(<alias>, "get_system"/"get_topology")(...)
            if (
                isinstance(f, ast.Call)
                and isinstance(f.func, ast.Name)
                and f.func.id == "getattr"
                and len(f.args) >= 2
                and isinstance(f.args[0], ast.Name)
                and f.args[0].id in aliases
                and isinstance(f.args[1], ast.Constant)
                and isinstance(f.args[1].value, str)
            ):
                if   f.args[1].value == "get_system":   self.found.add("system")
                elif f.args[1].value == "get_topology": self.found.add("topology")

            self.generic_visit(node)

    AliasCollector().visit(tree)
    depc = DepCollector()
    depc.visit(tree)
    return depc.found


__all__ = [
    "_discover_builtins_for_callable",
    "_ast_infer_builtins_for_callable",
]
