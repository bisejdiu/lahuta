# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     print(eval('"besian" + "sejdiu" + "@gmail.com"'))
#
"""Create a LahutaSystem directly from a model file using the fast model pathway."""

from pathlib import Path

from lahuta import InputType, LahutaSystem, logging
from lahuta.pipeline import InMemoryPolicy, Pipeline, PipelineContext
from lahuta.sources import FileSource


# fmt: off
def system_from_model_file(path: str | Path) -> LahutaSystem:
    sys = LahutaSystem(str(path), input_type=InputType.AlphaFold)
    ok = sys.build_topology()
    if not ok:
        raise RuntimeError("Failed to build topology from model file.")
    return sys


def pipeline_with_model_files(paths: list[str | Path]) -> None:
    p = Pipeline(FileSource([str(p) for p in paths]))
    p.params("system").is_model = True

    # Add a simple Python task that depends on the built-in 'system' node
    def inspect(ctx: PipelineContext):
        sys = ctx.get_system()
        return f"model_path={ctx.path}, n_atoms={sys.n_atoms}"

    p.add_task(name="inspect", task=inspect, depends=["system"], in_memory_policy=InMemoryPolicy.Keep)

    out = p.run(threads=1)
    for payload in out.raw("inspect"):
        logging.info(payload)


if __name__ == "__main__":
    logging.set_global_verbosity(logging.LogLevel.INFO)

    MODEL = Path(__file__).resolve().parents[3] / "core" / "data" / "fubi.cif"

    # Direct system creation from a model file
    sys = system_from_model_file(MODEL)
    top = sys.get_topology()
    logging.info(f"Direct model: atoms={sys.n_atoms}, residues={len(top.residues)}")

    # Pipeline over model files
    pipeline_with_model_files([MODEL])
