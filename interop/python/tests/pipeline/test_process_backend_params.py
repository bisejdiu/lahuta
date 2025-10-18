from __future__ import annotations

from pathlib import Path
from typing import Literal

import pytest

import lahuta as lxx
from lahuta.pipeline import InMemoryPolicy, Pipeline
from lahuta.sources import FileSource

# fmt: off
DATA_PATH = Path("core/data/models/AF-P0CL56-F1-model_v4.cif.gz")


def summarize_system(ctx):
    sys = ctx.get_system()
    top = ctx.get_topology()
    return {
        "is_model": bool(sys.is_model),
        "has_topology": bool(top is not None),
    }


def summarize_topology(ctx):
    topology = ctx.get_topology()
    return {
        "neighbors":    topology.is_computation_enabled(lxx.TopologyComputers.Neighbors),
        "bonds":        topology.is_computation_enabled(lxx.TopologyComputers.Bonds),
        "non_standard": topology.is_computation_enabled(lxx.TopologyComputers.NonStandardBonds),
        "residues":     topology.is_computation_enabled(lxx.TopologyComputers.Residues),
        "rings":        topology.is_computation_enabled(lxx.TopologyComputers.Rings),
        "atom_typing":  topology.is_computation_enabled(lxx.TopologyComputers.AtomTyping),
    }


@pytest.mark.parametrize("backend, extra", [("threads", {}), ("processes", {"processes": 2})],)
def test_process_backend_respects_params(backend: Literal["processes", "threads"], extra: dict[str, int]) -> None:
    if not DATA_PATH.exists():
        pytest.skip("ubi.cif fixture not available")

    source = FileSource([str(DATA_PATH)])
    pipeline = Pipeline(source)

    pipeline.params("system").is_model = True
    pipeline.params("topology").flags = lxx.TopologyComputers.Bonds | lxx.TopologyComputers.Residues
    pipeline.params("topology").atom_typing_method = lxx.AtomTypingMethod.GetContacts

    pipeline.add_task(
        name="system_info",
        task=summarize_system,
        depends=["topology"],
        in_memory_policy=InMemoryPolicy.Keep,
    )
    pipeline.add_task(
        name="topology_info",
        task=summarize_topology,
        depends=["topology"],
        in_memory_policy=InMemoryPolicy.Keep,
    )

    result = pipeline.run(threads=2, backend=backend, **extra)

    sys_info = result.json("system_info")[0]
    top_info = result.json("topology_info")[0]

    assert sys_info["is_model"] is True
    assert sys_info["has_topology"] is True

    assert top_info == {
        "neighbors":    False,
        "bonds":        True,
        "non_standard": False,
        "residues":     True,
        "rings":        False,
        "atom_typing":  False,
    }
