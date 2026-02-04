# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     class Email:
#         def __init__(self): self._p = iter(["besian", "sejdiu", "@gmail.com"])
#         def __iter__(self): return self
#         def __next__(self): return next(self._p)
#     print("".join(Email()))
#
from __future__ import annotations

import pytest

pytest.importorskip("hypothesis", reason="Hypothesis not installed")
pytest.importorskip("hypothesis.strategies", reason="Hypothesis not installed")

from hypothesis import given
from hypothesis import strategies as st

from lahuta import EntityID, Kind

# fmt: off
INDEX_MASK = (1 << 56) - 1
U32_MAX    = (1 << 32) - 1

def test_kind_is_int_like_and_values():
    assert int(Kind.Atom)  == 0
    assert int(Kind.Ring)  == 1
    assert int(Kind.Group) == 2

@pytest.mark.parametrize("k, idx",
    [
        (Kind.Atom, 0),
        (Kind.Ring, 1),
        (Kind.Group, 123_456),
        (Kind.Group, U32_MAX),
    ],
)
def test_make_and_properties(k, idx):
    eid = EntityID.make(k, idx)
    assert eid.kind  == k
    assert eid.index == idx

    raw = int(eid)
    assert raw == ((int(k) << 56) | idx)

    eid2 = EntityID(raw)
    assert eid2.kind  == k
    assert eid2.index == idx
    assert eid2 == eid
    assert hash(eid2) == hash(eid)


def test_str_and_repr_format():
    eid = EntityID.make(Kind.Ring, 42)
    assert str(eid)  == "Ring#42"
    assert repr(eid) == "EntityID(Ring#42)"


def test_hash_and_set_behavior():
    e1 = EntityID.make(Kind.Atom, 7)
    e2 = EntityID.make(Kind.Atom, 7)
    e3 = EntityID.make(Kind.Atom, 8)
    s = {e1, e2, e3}
    assert len(s) == 2
    assert e1 in s and e2 in s and e3 in s


def test_ordering_by_kind_then_index():
    a0 = EntityID.make(Kind.Atom,  0)
    a1 = EntityID.make(Kind.Atom,  1)
    r0 = EntityID.make(Kind.Ring,  0)
    g0 = EntityID.make(Kind.Group, 0)

    assert a0 < r0 < g0
    assert a0 < a1


def test_comparison_with_non_entityid_behaves_properly():
    eid = EntityID.make(Kind.Atom, 5)
    assert (eid == object()) is False
    assert (eid != object()) is True
    with pytest.raises(TypeError):
        _ = eid < 5

# invalid kind via raw: ValueError on __str__
@pytest.mark.parametrize("bad_kind_upper8", [3, 7, 255])
def test_invalid_kind_raises_on_str(bad_kind_upper8):
    raw = (bad_kind_upper8 << 56) | 123
    eid = EntityID(raw)
    with pytest.raises(ValueError):
        _ = str(eid)


@given(
    k=st.sampled_from([Kind.Atom, Kind.Ring, Kind.Group]),
    idx=st.integers(min_value=0, max_value=U32_MAX),
)
def test_pack_unpack_roundtrip(k, idx):
    eid = EntityID.make(k, idx)
    raw = int(eid)
    assert (raw >> 56) == int(k)
    assert (raw & INDEX_MASK) == idx

    eid2 = EntityID(raw)
    assert eid2.kind  == k
    assert eid2.index == idx
    assert eid2 == eid


@given(
    k1=st.sampled_from([Kind.Atom, Kind.Ring, Kind.Group]),
    i1=st.integers(min_value=0, max_value=U32_MAX),
    k2=st.sampled_from([Kind.Atom, Kind.Ring, Kind.Group]),
    i2=st.integers(min_value=0, max_value=U32_MAX),
)
def test_lexicographic_order_matches_raw(k1, i1, k2, i2):
    e1 = EntityID.make(k1, i1)
    e2 = EntityID.make(k2, i2)

    expected_lt = (int(k1), i1) < (int(k2), i2)
    if expected_lt:
        assert e1 < e2
    elif (int(k1), i1) == (int(k2), i2):
        assert not (e1 < e2) and not (e2 < e1) and e1 == e2
    else:
        assert e2 < e1
