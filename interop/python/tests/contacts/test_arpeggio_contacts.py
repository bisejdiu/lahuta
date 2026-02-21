# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     def f(**kw):
#         return kw["f"] + kw["l"] + "@" + kw["d"]
#     print(f(**{"f": "besian", "l": "sejdiu", "d": "gmail.com"}))
#
from __future__ import annotations

import os
import re
from collections import Counter
from pathlib import Path

import pytest

from lahuta import AtomTypingMethod, ContactProvider
from lahuta.db import LahutaDB
from lahuta.pipeline import Pipeline
from lahuta.pipeline.tasks import ContactTask
from lahuta.sources import DatabaseHandleSource, FileSource


def type_counts(rec: dict) -> Counter:
    return Counter(str(c.get("type", "")) for c in rec.get("contacts", []))


def pair_multiset(rec: dict) -> Counter:
    pairs = (
        "|".join(
            sorted(
                [
                    _canon_endpoint(str(c.get("lhs", ""))),
                    _canon_endpoint(str(c.get("rhs", ""))),
                ]
            )
        )
        for c in rec.get("contacts", [])
    )
    return Counter(pairs)


# Canonicalize group repr
_GROUP_RE = re.compile(
    r"^\(\s*(?P<atoms>.+?)\s*\)\s*-\s*(?P<resseq>[^-\s]+)\s*-\s*(?P<resname>[^-\s]+)\s*-\s*(?P<chain>\S+)\s*$"
)


def _canon_endpoint(s: str) -> str:
    s = " ".join(s.split())
    m = _GROUP_RE.match(s)
    if not m:
        return s
    atoms = [tok.strip() for tok in m.group("atoms").split(",") if tok.strip()]
    atoms.sort()
    return f"({', '.join(atoms)})-{m.group('resseq')}-{m.group('resname')}-{m.group('chain')}"


@pytest.mark.parametrize(
    "model_basename",
    [
        "AF-P0CL56-F1-model_v4.cif.gz",
        "AF-Q57552-F1-model_v4.cif.gz",
    ],
)
def test_arpeggio_contacts_file_vs_db(tmp_path: Path, model_basename: str) -> None:
    data_dir = Path("core/data/models")
    assert data_dir.exists(), "Models directory missing"

    model_path = data_dir / model_basename
    assert model_path.exists(), f"Missing model file: {model_path}"

    p_file = Pipeline(FileSource(str(model_path)))
    p_file.add_task(name="contacts_file", task=ContactTask(provider=ContactProvider.Arpeggio))
    out_file = p_file.run(threads=1)
    recs_file = out_file.to_dict("contacts_file")
    assert isinstance(recs_file, list) and len(recs_file) == 1
    rec_file = recs_file[0]

    db_path = tmp_path / "models_lmdb"
    ldb = LahutaDB.create_from_directory(data_dir, db_path, ext=".cif.gz", recursive=False, batch=50, threads=1)
    p_db = Pipeline(DatabaseHandleSource(ldb._db, batch=64))

    # Ensure Arpeggio atom typing is selected before building topology in model mode
    p_db.params("topology").atom_typing_method = AtomTypingMethod.Arpeggio
    p_db.add_task(name="contacts_db", task=ContactTask(provider=ContactProvider.Arpeggio))
    out_db = p_db.run(threads=1)
    recs_db = out_db.to_dict("contacts_db")
    assert isinstance(recs_db, list) and len(recs_db) >= 1

    rec_db = next((r for r in recs_db if os.path.basename(str(r.get("file_path", ""))) == model_basename), None)
    assert rec_db is not None, f"DB output missing record for {model_basename}"

    assert str(rec_file.get("provider", "")).lower() == "arpeggio"
    assert str(rec_db.get("provider", "")).lower() == "arpeggio"
    assert isinstance(rec_file.get("contacts", []), list)
    assert isinstance(rec_db.get("contacts", []), list)

    tc_file = type_counts(rec_file)
    tc_db = type_counts(rec_db)
    assert tc_file == tc_db, f"Contact-type counts differ for {model_basename}: file={tc_file} db={tc_db}"

    pm_file = pair_multiset(rec_file)
    pm_db = pair_multiset(rec_db)
    assert pm_file == pm_db, (
        f"Participant pairs differ for {model_basename}: file_only={len(pm_file - pm_db)} db_only={len(pm_db - pm_file)}"
    )
