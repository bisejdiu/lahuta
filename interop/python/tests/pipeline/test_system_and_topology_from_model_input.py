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

from pathlib import Path


def test_direct_model_then_pipeline_fast_path(data_dir, run_child) -> None:
    """
    Both LahutaSystem and Pipeline will use the same thread and thus
    the same thread-local cache.
    """
    p = Path(data_dir) / "fubi.cif"
    model = str(p)
    child = f"""
from lahuta import LahutaSystem
from lahuta.pipeline import Pipeline
from lahuta.sources import FileSource

# Direct fast-path system
sys = LahutaSystem.from_model_file({model!r})
ok = sys.build_topology()
assert ok is True
assert sys.has_topology_built is True or sys.has_topology_built()

# Pipeline fast-path in same process
p = Pipeline(FileSource({model!r}))
p.params("system").is_model = True

def inspect(ctx) -> str: # our type inference cannot infer str in this context
    s = ctx.get_system()
    topo = ctx.get_topology()
    assert topo is not None
    payload = ctx.model_payload
    assert payload is not None
    dssp = payload.dssp
    residue_bf = payload.bfactors
    ctx_bf = ctx.bfactors
    assert isinstance(dssp, list) and isinstance(residue_bf, list)
    assert len(dssp) == len(residue_bf) and len(dssp) > 0
    assert ctx_bf is None or isinstance(ctx_bf, list)
    meta = payload.metadata
    assert meta is None or isinstance(meta, dict)
    return f"ok {{ctx.path}} {{int(s.n_atoms)}}"

p.add_task(name="inspect", task=inspect, depends=["topology"])
out = p.run(threads=1)
print("OK", out.get("inspect", []))
"""

    res = run_child(child)
    if res.returncode != 0:
        raise AssertionError(
            f"child failed (returncode={res.returncode})\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        )
    assert "OK" in res.stdout
