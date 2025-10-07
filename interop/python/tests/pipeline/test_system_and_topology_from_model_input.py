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
import lahuta as lxx
from lahuta.pipeline import Pipeline
from lahuta.sources import FileSource

# Direct fast-path system
sys = lxx.LahutaSystem.from_model_file({model!r})
ok = sys.build_topology()
assert ok is True
assert sys.has_topology_built is True or sys.has_topology_built()

# Pipeline fast-path in same process
p = Pipeline(FileSource({model!r}))
p.params("system").is_model = True

def inspect(ctx):
    s = ctx.get_system()
    return f"ok {{ctx.path}} {{int(s.n_atoms)}}"

p.add_task(name="inspect", task=inspect, depends=["system"])
out = p.run(threads=1)
print("OK", out.get("inspect", []))
"""

    res = run_child(child)
    if res.returncode != 0:
        raise AssertionError(
            f"child failed (returncode={res.returncode})\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        )
    assert "OK" in res.stdout
