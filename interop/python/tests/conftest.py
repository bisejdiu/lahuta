from pathlib import Path
from typing import Callable

import pytest


def _resolve_data_file(filename: str) -> Path:
    here = Path(__file__).resolve()
    p = here
    for _ in range(12):
        candidate = p.parent / "data" / filename if p.name == "tests" else p / "data" / filename
        if candidate.exists():
            return candidate
        direct = p / "data" / filename
        if direct.exists():
            return direct
        p = p.parent
    alt = Path.cwd() / "data" / filename
    if alt.exists():
        return alt
    pytest.skip(f"Missing data/{filename}")


@pytest.fixture(scope="session")
def data_dir() -> Path:
    here = Path(__file__).resolve()
    for p in [here.parent, *here.parents]:
        if (p / "core").is_dir() and (p / "interop" / "python").is_dir():
            root = (
                p.parent
                if p.name == "interop" and (p.parent / "core").is_dir() and (p.parent / "interop" / "python").is_dir()
                else p
            )
            d = root / "core" / "data"
            if d.is_dir():
                return d
    cand = Path.cwd() / "core" / "data"
    if cand.is_dir():
        return cand
    pytest.skip("Missing core/data")


@pytest.fixture
def test_files(data_dir: Path) -> list[str]:
    files: list[str] = []
    for pat in ("*.cif", "*.cif.gz"):
        files.extend([str(p) for p in sorted(data_dir.glob(pat))])
    if not files:
        pytest.skip(f"No CIF files in {data_dir}")
    return files[:3] if len(files) >= 3 else files


@pytest.fixture(scope="session")
def data_path() -> Callable[[str], Path]:
    def _inner(name: str) -> Path:
        return _resolve_data_file(name)

    return _inner


@pytest.fixture(scope="session")
def ubi_cif(data_path: Callable[[str], Path]) -> Path:
    return data_path("ubi.cif")


@pytest.fixture(scope="session")
def luni(ubi_cif: Path):
    import lahuta as lxx

    return lxx.LahutaSystem(str(ubi_cif))


@pytest.fixture(scope="session")
def run_child() -> Callable[[str], object]:
    # Execute a Python code snippet in a subprocess and capture output
    import os
    import subprocess
    import sys
    from pathlib import Path as _Path

    def _inner(code: str) -> subprocess.CompletedProcess:
        return subprocess.run(
            [sys.executable, "-c", code],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env={**os.environ},
            cwd=str(_Path.cwd()),
        )

    return _inner
