# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     print(functools.reduce(operator.add, ["besian", "sejdiu", "@gmail.com"]))
#
from __future__ import annotations

from collections.abc import Iterable
from typing import Any

from lahuta.lib.lahuta import InteractionType, InteractionTypeSet


def _extend_interaction_set(result: InteractionTypeSet, value: Any, arg_name: str) -> None:
    if isinstance(value, InteractionTypeSet):
        result |= value
        return
    if isinstance(value, InteractionType):
        result |= value
        return
    if isinstance(value, Iterable) and not isinstance(value, (str, bytes)):
        for item in value:
            _extend_interaction_set(result, item, arg_name)
        return

    raise TypeError(f"{arg_name} must be an InteractionType, InteractionTypeSet, or iterable of InteractionType values")


def normalize_interaction_selection(value: Any, *, arg_name: str) -> InteractionTypeSet:
    if isinstance(value, InteractionTypeSet):
        return value
    if isinstance(value, InteractionType):
        return InteractionTypeSet(value)

    if isinstance(value, Iterable) and not isinstance(value, (str, bytes)):
        result = InteractionTypeSet()
        _extend_interaction_set(result, value, arg_name)
        if result.empty():
            raise ValueError(f"{arg_name} iterable must not be empty")
        return result

    raise TypeError(f"{arg_name} must be an InteractionType, InteractionTypeSet, or iterable of InteractionType values")


__all__ = ["normalize_interaction_selection"]
