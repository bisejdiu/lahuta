from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Generic, TypeVar

from ..lib import lahuta as lxx
from ..lib.lahuta import AtomRec, GroupRec, RingRec

TRec_co = TypeVar("TRec_co", AtomRec, RingRec, GroupRec, covariant=True)


# fmt: off
@dataclass(frozen=True)
class Selector(Generic[TRec_co]):
    """Selection for a given entity kind with an optional predicate."""

    kind: lxx.Kind
    predicate: Callable[[TRec_co], bool] | None = None


def atoms(predicate:  Callable[[AtomRec], bool] | None = None) ->  Selector[AtomRec]:
    return Selector(kind=lxx.Kind.Atom, predicate=predicate)


def rings(predicate:  Callable[[RingRec], bool] | None = None) ->  Selector[RingRec]:
    return Selector(kind=lxx.Kind.Ring, predicate=predicate)


def groups(predicate: Callable[[GroupRec], bool] | None = None) -> Selector[GroupRec]:
    return Selector(kind=lxx.Kind.Group, predicate=predicate)


__all__ = [
    "Selector",
    "atoms",
    "rings",
    "groups",
]
