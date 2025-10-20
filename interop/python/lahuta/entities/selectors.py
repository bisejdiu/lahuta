from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Generic, TypeVar

from ..lib import lahuta as lxx
from ..lib.lahuta import AtomRec, GroupRec, RingRec

TRec_co = TypeVar("TRec_co", AtomRec, RingRec, GroupRec, covariant=True)


# fmt: off
@dataclass(frozen=True)
class Selector(Generic[TRec_co]):
    """Selection for a given entity kind with an optional selector."""

    kind: lxx.Kind
    selector: Callable[[TRec_co], bool] | None = None


def atoms(selector:  Callable[[AtomRec], bool] | None = None) ->  Selector[AtomRec]:
    return Selector(kind=lxx.Kind.Atom, selector=selector)


def rings(selector:  Callable[[RingRec], bool] | None = None) ->  Selector[RingRec]:
    return Selector(kind=lxx.Kind.Ring, selector=selector)


def groups(selector: Callable[[GroupRec], bool] | None = None) -> Selector[GroupRec]:
    return Selector(kind=lxx.Kind.Group, selector=selector)


__all__ = [
    "Selector",
    "atoms",
    "rings",
    "groups",
]
