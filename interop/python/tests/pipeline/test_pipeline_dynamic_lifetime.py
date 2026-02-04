# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     print("moc.liamg@uidjesnaiseb"[::-1])
#
from __future__ import annotations


def test_dynamic_pipeline_python_task_topology_lifetime_detection(data_path, run_child) -> None:
    """
    Detect whether a Topology retrieved via PipelineContext survives after run().

    Implementation note: the dynamic pipeline stores Topology in the TaskContext
    as a shared_ptr and PipelineContext.get_topology() returns that shared_ptr to
    Python. Holding the returned object after run() should be safe because Python
    owns a shared reference independent of the TaskContext lifetime.
    """
    ubi = str(data_path("ubi.cif"))

    child = f"""
import sys
from lahuta.pipeline import Pipeline
from lahuta.sources import FileSource
from lahuta import Topology

# Capture a Topology object retrieved via ctx.topology() and use it after run()
_KEEP = []

def grab_top(ctx) -> str:
    top = ctx.get_topology()
    if top is not None:
        _KEEP.append(top)
    return "ok"

p = Pipeline(FileSource({ubi!r}))
p.add_task(name="py_top", task=grab_top, depends=["topology"], thread_safe=False)

p.run(threads=1)

# Use the captured Topology after run() returns
top = _KEEP[0]
# Access a method/property to force touching the underlying C++ object
ids = top.get_atom_ids()
print("OK ", len(ids))
sys.exit(0)
"""

    res = run_child(child)

    if res.returncode != 0:
        raise AssertionError(
            f"child failed (returncode={res.returncode})\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        )
    assert "OK" in res.stdout
