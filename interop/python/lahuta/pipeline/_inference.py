# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     class Email:
#         def __init__(self, s=""): self.s = s
#         def __matmul__(self, other): return Email(self.s + other)
#         def __str__(self): return self.s
#     print(str(Email() @ "besian" @ "sejdiu" @ "@gmail.com"))
#
from __future__ import annotations

import ast
import inspect
from typing import Any, Callable, get_args, get_origin, get_type_hints, Literal

from .types import OutputFormat

_JSON_SCALARS = (type(None), bool, int, float)  # str means TEXT

# We want to detect:
# - TEXT (strong): str, Literal[str,...], Union[str,None], Optional[str]
# - JSON (simple): dict, list, Mapping, Sequence, Union[non-str scalars], Optional[non-str scalar]
# - Fallback: JSON


def _unwrap_annotated(tp: Any) -> Any:
    origin = get_origin(tp)
    if origin is None:
        return tp
    try:
        from typing import Annotated
    except Exception:
        Annotated = None
    if Annotated is not None and origin is Annotated:
        args = get_args(tp)
        return args[0] if args else tp
    return tp


def _is_none_type(tp: Any) -> bool:
    return tp is type(None)


def _is_json_container_origin(origin: Any) -> bool:
    if origin in (list, dict):
        return True
    try:
        from collections.abc import Mapping, MutableMapping, Sequence, MutableSequence

        return origin in (Mapping, MutableMapping, Sequence, MutableSequence)
    except Exception:
        return False


def _classify_annotation(tp: Any) -> OutputFormat | None:
    """Strong TEXT detection, simple JSON detection."""
    if tp is None:
        return None

    tp = _unwrap_annotated(tp)
    origin = get_origin(tp)
    args = get_args(tp)

    # Direct primitives
    if tp is str:
        return OutputFormat.TEXT
    if tp is bytes:
        return OutputFormat.BINARY
    if tp in _JSON_SCALARS and tp is not str:
        return OutputFormat.JSON

    # Literal[...] (all strings -> TEXT, all non-strings -> JSON)
    try:
        if origin is Literal and args:
            all_str = all(isinstance(a, str) for a in args)
            all_non_str = all(not isinstance(a, str) for a in args)
            if all_str:
                return OutputFormat.TEXT
            if all_non_str:
                return OutputFormat.JSON
            return None
    except Exception:
        pass

    # Containers / generics
    if tp in (dict, list) or _is_json_container_origin(origin):
        return OutputFormat.JSON

    # Union/Optional
    if origin is getattr(__import__("typing"), "Union", None):
        arm_kinds: set[OutputFormat | None] = set()
        for a in args:
            a = _unwrap_annotated(a)
            if a is str:
                arm_kinds.add(OutputFormat.TEXT)
            elif _is_none_type(a):
                arm_kinds.add(None)
            else:
                arm_kinds.add(_classify_annotation(a))
        # TEXT | None -> TEXT
        if arm_kinds in ({OutputFormat.TEXT, None}, {OutputFormat.TEXT}, {None}):
            return OutputFormat.TEXT
        # All JSON-ish (maybe with None) -> JSON
        ak = {k for k in arm_kinds if k is not None}
        if ak == {OutputFormat.JSON}:
            return OutputFormat.JSON
        return None

    return None


def _infer_from_annotation(task: Callable[..., Any]) -> OutputFormat | None:
    try:
        hints = get_type_hints(task, include_extras=True)
    except Exception:
        return None
    rt = hints.get("return")
    return _classify_annotation(rt) if rt is not None else None


def _resolve_name_value(task: Callable[..., Any], name: str) -> Any:
    """
    Resolve Name to a bound runtime constant via closure/globals (no execution).
    Used only for TEXT checks.
    """
    try:
        cv = inspect.getclosurevars(task)
    except Exception:
        cv = None
    if cv:
        if name in cv.nonlocals:
            return cv.nonlocals[name]
        if name in cv.globals:
            return cv.globals[name]
        if name in cv.builtins:
            return cv.builtins[name]
    try:
        g = getattr(task, "__globals__", None) or {}
        if name in g:
            return g[name]
    except Exception:
        pass
    return inspect._empty


def _is_textlike_expr(node: ast.AST, *, task: Callable[..., Any]) -> bool:
    # "literal"
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return True
    # f"..."
    if isinstance(node, ast.JoinedStr):
        return True
    # "a" + x
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        return _is_textlike_expr(node.left, task=task) or _is_textlike_expr(node.right, task=task)
    # "a" % x
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mod):
        return _is_textlike_expr(node.left, task=task)
    # "a".format(...), "".join(...)
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
        if node.func.attr in {"format", "join"}:
            return _is_textlike_expr(node.func.value, task=task)
    # str(...)
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "str":
        return True
    # Name resolved via closure/globals
    if isinstance(node, ast.Name):
        val = _resolve_name_value(task, node.id)
        return isinstance(val, str)
    # IfExp/BoolOp: any arm textlike?
    if isinstance(node, ast.IfExp):
        return _is_textlike_expr(node.body, task=task) or _is_textlike_expr(node.orelse, task=task)
    if isinstance(node, ast.BoolOp):
        return any(_is_textlike_expr(v, task=task) for v in node.values)
    return False


def _is_jsonlike_expr_simple(node: ast.AST) -> bool:
    """
    Minimal JSON shape detection: non-string constants, dict/list literals. dict()/list() calls.
    No name resolution, no tuples/sets, intentionally conservative.
    """
    if isinstance(node, ast.Constant) and not isinstance(node.value, str):
        return True
    if isinstance(node, (ast.Dict, ast.List)):
        return True
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in {"dict", "list"}:
        return True
    return False


def _find_target_func_node(mod: ast.Module, target_name: str) -> ast.AST | None:
    for n in ast.walk(mod):
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name == target_name:
            return n
    lambdas = [n for n in ast.walk(mod) if isinstance(n, ast.Lambda)]
    return lambdas[0] if len(lambdas) == 1 else None


def _infer_from_ast(task: Callable[..., Any]) -> OutputFormat | None:
    try:
        src = inspect.getsource(task)
    except (OSError, TypeError):
        return None

    import textwrap

    src = textwrap.dedent(src)
    try:
        mod = ast.parse(src)
    except SyntaxError:
        return None

    node = _find_target_func_node(mod, getattr(task, "__name__", ""))
    if node is None:
        return None

    saw_text = False
    saw_json = False
    saw_only_none = True

    def _check(val: ast.expr | None) -> None:
        nonlocal saw_text, saw_json, saw_only_none
        if val is None or (isinstance(val, ast.Constant) and val.value is None):
            return
        saw_only_none = False
        if _is_textlike_expr(val, task=task):
            saw_text = True
            return
        if _is_jsonlike_expr_simple(val):
            saw_json = True

    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
        for n in ast.walk(node):
            if isinstance(n, ast.Return):
                _check(n.value)
    elif isinstance(node, ast.Lambda):
        _check(node.body)

    if saw_text and not saw_json:
        return OutputFormat.TEXT
    if saw_json and not saw_text:
        return OutputFormat.JSON
    if saw_only_none:
        return OutputFormat.JSON
    return None


def infer_python_task_format(task: Callable[..., Any]) -> OutputFormat:
    """
    - Strong TEXT via hints/AST, with Name-resolution for strings
    - Simple JSON via hints/AST, no exhaustive edge cases
    - Default is JSON
    """
    ann = _infer_from_annotation(task)
    if ann is not None:
        return ann

    ast_guess = _infer_from_ast(task)
    if ast_guess is not None:
        return ast_guess

    try:
        import logging

        logging.info(
            "Could not infer return type for task %r. Defaulting to JSON. "
            "Add type annotation (-> str or -> dict) for explicit control.",
            getattr(task, "__name__", task),
        )
    except Exception:
        pass

    return OutputFormat.JSON
