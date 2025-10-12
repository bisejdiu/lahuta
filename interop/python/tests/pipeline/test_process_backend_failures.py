from __future__ import annotations

import os
import tempfile
import textwrap

import pytest

from lahuta.pipeline import InMemoryPolicy, Pipeline


def _run_child(code: str) -> tuple[int, str, str]:
    import subprocess
    import sys

    dedented = textwrap.dedent(code).strip() + "\n"
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write(dedented)
        script = fh.name
    try:
        proc = subprocess.run(
            [sys.executable, script],
            capture_output=True,
            text=True,
            env=dict(os.environ),
            cwd=os.getcwd(),
        )
        return proc.returncode, proc.stdout, proc.stderr
    finally:
        try:
            os.remove(script)
        except OSError:
            pass


def test_process_backend_timeout_emits_error(data_path) -> None:
    data_file = data_path("models/AF-P0CL56-F1-model_v4.cif.gz")

    child = f"""
import json
import time
from multiprocessing import freeze_support

from lahuta.pipeline import Pipeline, InMemoryPolicy
from lahuta.sources import FileSource


def slow_task(_):
    time.sleep(1.0)
    return "slow"


def main():
    p = Pipeline(FileSource([{str(data_file)!r}]))
    p.add_task(name="slow", task=slow_task, in_memory_policy=InMemoryPolicy.Keep)

    res = p.run(threads=2, backend="processes", processes=1, process_timeout=0.1)
    payload = res.raw("slow")[0]
    data = json.loads(payload)
    assert data["error"]["message"] == "Timeout", payload
    print("OK")


if __name__ == "__main__":
    freeze_support()
    main()
"""

    rc, out, err = _run_child(child)
    if rc != 0:
        raise AssertionError(f"child failed (returncode={rc})\nSTDOUT:\n{out}\nSTDERR:\n{err}")
    assert "OK" in out


def test_process_backend_worker_crash_emits_error(data_path) -> None:
    data_file = data_path("models/AF-P0CL56-F1-model_v4.cif.gz")

    child = f"""
import json
import os
from multiprocessing import freeze_support
from lahuta.pipeline import Pipeline, InMemoryPolicy
from lahuta.sources import FileSource


def crash_task(_):
    os._exit(1)


def main():
    p = Pipeline(FileSource([{str(data_file)!r}]))
    p.add_task(name="crash", task=crash_task, in_memory_policy=InMemoryPolicy.Keep, depends=[])

    # Note: Python's multiprocessing doesn't distinguish worker crashes from timeouts well.
    # When a worker calls os._exit(), the pool.apply_async().get(timeout) will wait the
    # full timeout before raising TimeoutError. Use a short timeout to make the test fast.
    res = p.run(threads=2, backend="processes", processes=1, process_timeout=0.3)
    payload = res.raw("crash")[0]
    data = json.loads(payload)
    # Worker crash manifests as timeout in Python's multiprocessing
    assert "error" in data, payload

    # double braces are necessary for proper f-string formatting
    assert data["error"]["message"] in ["Timeout", "Worker crashed"], f'Got: {{data["error"]["message"]}}'
    print("OK")


if __name__ == "__main__":
    freeze_support()
    main()
"""

    rc, out, err = _run_child(child)
    if rc != 0:
        raise AssertionError(f"child failed (returncode={rc})\nSTDOUT:\n{out}\nSTDERR:\n{err}")
    assert "OK" in out


def test_process_backend_worker_crash_allows_subsequent_runs(data_path) -> None:
    data_file = data_path("models/AF-P0CL56-F1-model_v4.cif.gz")

    child = f"""
import json
import os
from multiprocessing import freeze_support
from lahuta.pipeline import Pipeline, InMemoryPolicy
from lahuta.sources import FileSource


def crash_task(_):
    os._exit(1)


def ok_task(_):
    return "ok"


def main():
    crash_pipeline = Pipeline(FileSource([{str(data_file)!r}]))
    crash_pipeline.add_task(name="crash", task=crash_task, in_memory_policy=InMemoryPolicy.Keep, depends=[])

    crash_result = crash_pipeline.run(threads=2, backend="processes", processes=1, process_timeout=0.3)
    crash_payload = crash_result.raw("crash")[0]
    crash_data = json.loads(crash_payload)
    assert "error" in crash_data, crash_payload
    assert crash_data["error"]["message"] in ["Timeout", "Worker crashed"], crash_payload

    ok_pipeline = Pipeline(FileSource([{str(data_file)!r}]))
    ok_pipeline.add_task(name="ok", task=ok_task, in_memory_policy=InMemoryPolicy.Keep, depends=[])

    ok_result = ok_pipeline.run(threads=2, backend="processes", processes=1, process_timeout=5.0)
    ok_payload = ok_result.raw("ok")[0]
    assert ok_payload == "ok", ok_payload
    print("OK")


if __name__ == "__main__":
    freeze_support()
    main()
"""

    rc, out, err = _run_child(child)
    if rc != 0:
        raise AssertionError(f"child failed (returncode={rc})\nSTDOUT:\n{out}\nSTDERR:\n{err}")
    assert "OK" in out


def test_process_backend_partial_worker_crash_allows_subsequent_runs(data_path) -> None:
    data_file = data_path("models/AF-P0CL56-F1-model_v4.cif.gz")

    child = f"""
import json
import os
from multiprocessing import freeze_support
from lahuta.pipeline import Pipeline, InMemoryPolicy
from lahuta.sources import FileSource


def crash_task(_):
    os._exit(1)


def ok_task(_):
    return "ok"


def main():
    crash_pipeline = Pipeline(FileSource([{str(data_file)!r}]))
    crash_pipeline.add_task(name="crash", task=crash_task, in_memory_policy=InMemoryPolicy.Keep, depends=[])

    crash_result = crash_pipeline.run(threads=2, backend="processes", processes=2, process_timeout=0.3)
    crash_payload = crash_result.raw("crash")[0]
    crash_data = json.loads(crash_payload)
    assert "error" in crash_data, crash_payload
    assert crash_data["error"]["message"] in ["Timeout", "Worker crashed"], crash_payload

    ok_pipeline = Pipeline(FileSource([{str(data_file)!r}]))
    ok_pipeline.add_task(name="ok", task=ok_task, in_memory_policy=InMemoryPolicy.Keep, depends=[])

    ok_result = ok_pipeline.run(threads=2, backend="processes", processes=2, process_timeout=5.0)
    ok_payload = ok_result.raw("ok")[0]
    assert ok_payload == "ok", ok_payload
    print("OK")


if __name__ == "__main__":
    freeze_support()
    main()
"""

    rc, out, err = _run_child(child)
    if rc != 0:
        raise AssertionError(f"child failed (returncode={rc})\nSTDOUT:\n{out}\nSTDERR:\n{err}")
    assert "OK" in out


def test_process_backend_disallows_md_sources(data_path) -> None:
    md_structure = data_path("simulationdatabase/lysozyme.gro")
    md_xtc = data_path("simulationdatabase/lysozyme.xtc")

    from lahuta.lib import lahuta as _lib

    md_source = _lib.pipeline.sources.MdTrajectoriesSource(
        [
            {
                "structure": str(md_structure),
                "xtcs": [str(md_xtc)],
                "id": "traj",
            }
        ]
    )

    pipeline = Pipeline(md_source)  # type: ignore
    pipeline.add_task(
        name="noop",
        task=lambda _: "ok",
        in_memory_policy=InMemoryPolicy.Keep,
    )

    with pytest.raises(RuntimeError, match="backend='processes' is not supported"):
        pipeline.run(backend="processes", processes=1)


def test_process_backend_executes_notebook_defined_task_in_worker(data_path) -> None:
    """Test that notebook-defined tasks (no importable module) can be executed via process backend."""
    data_file = data_path("models/AF-P0CL56-F1-model_v4.cif.gz")

    child = f"""
import os
from multiprocessing import freeze_support
from lahuta.pipeline import Pipeline, InMemoryPolicy
from lahuta.sources import FileSource


def main():
    parent_pid = os.getpid()

    # Simulate notebook task: defined in __main__ with no importable module/qualname
    def notebook_task(ctx):
        import os
        return {{"worker_pid": os.getpid()}}

    pipeline = Pipeline(FileSource([{str(data_file)!r}]))
    pipeline.add_task(name="pid", task=notebook_task, in_memory_policy=InMemoryPolicy.Keep, depends=[])

    result = pipeline.run(threads=2, backend="processes", processes=1, process_timeout=10.0)
    payload = result.json("pid")[0]
    worker_pid = payload.get("worker_pid")
    assert worker_pid is not None, f"Expected worker_pid in payload, got: {{payload}}"
    assert worker_pid != parent_pid, f"Expected worker PID != parent (got {{worker_pid}}, parent={{parent_pid}})"
    print("OK")


if __name__ == "__main__":
    freeze_support()
    main()
"""

    rc, out, err = _run_child(child)
    if rc != 0:
        raise AssertionError(f"child failed (returncode={rc})\nSTDOUT:\n{out}\nSTDERR:\n{err}")
    assert "OK" in out


def test_process_backend_notebook_task_with_topology(data_path) -> None:
    """
    Test notebook-defined tasks can access topology and return structured data.

    This complements test_process_backend_executes_notebook_defined_task_in_worker by:
    - Testing topology access (not just PID)
    - Validating actual data processing
    - Running directly in test process (no subprocess isolation)

    The task is defined in the test module scope (not __main__), simulating
    a notebook cell definition that cannot be imported elsewhere.
    """
    data_file = data_path("models/AF-P0CL56-F1-model_v4.cif.gz")

    from lahuta.sources import FileSource

    # Define task in test scope - simulates notebook cell definition
    # This cannot be imported by module/qualname from other modules
    def notebook_analysis_task(ctx):
        """Task defined in test scope - must be serialized via cloudpickle."""
        import os

        worker_pid = os.getpid()
        topology = ctx.get_topology()

        if topology is None:
            return {"error": "Failed to build topology", "worker_pid": worker_pid}

        atom_ids = list(topology.get_atom_ids())
        return {
            "worker_pid": worker_pid,
            "path": ctx.path,
            "n_atoms": len(atom_ids),
            "has_topology": True,
        }

    # Verify task is defined in test module (like notebook __main__)
    assert notebook_analysis_task.__module__ == "test_process_backend_failures", (
        f"Expected task in test module, got {notebook_analysis_task.__module__}"
    )

    pipeline = Pipeline(FileSource([str(data_file)]))
    pipeline.add_task(
        name="analysis",
        task=notebook_analysis_task,
        in_memory_policy=InMemoryPolicy.Keep,
    )

    # If cloudpickle not available or serialization fails, this should fall back to
    # threaded backend with a warning
    result = pipeline.run(threads=2, backend="processes", processes=2, process_timeout=10.0)
    payload = result.json("analysis")[0]

    # Validate execution in worker process (may be parent if fallback occurred)
    worker_pid = payload.get("worker_pid")
    assert worker_pid is not None, f"Expected worker_pid in payload, got: {payload}"

    # Validate topology access worked
    assert payload.get("has_topology") is True, "Expected topology to be built"
    assert payload.get("n_atoms", 0) > 0, f"Expected positive atom count, got: {payload.get('n_atoms')}"
    assert "error" not in payload, f"Task failed with error: {payload.get('error')}"

    # Validate data integrity
    assert isinstance(payload.get("path"), str), "Expected path to be string"
    assert data_file.name in payload["path"], f"Expected {data_file.name} in path {payload['path']}"
