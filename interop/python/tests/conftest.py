import os
from pathlib import Path
from typing import Callable

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--lahuta-error",
        action="store_true",
        default=False,
        help="Set Lahuta logging to ERROR level",
    )
    parser.addoption(
        "--lahuta-info",
        action="store_true",
        default=False,
        help="Set Lahuta logging to INFO level",
    )
    parser.addoption(
        "--lahuta-debug",
        action="store_true",
        default=False,
        help="Set Lahuta logging to DEBUG level",
    )


def pytest_configure(config):
    from lahuta import logging
    from lahuta.logging import LogLevel

    # debug > info > error
    if config.option.lahuta_debug:
        logging.set_global_verbosity(LogLevel.DEBUG)
        print("[pytest] Lahuta logging: DEBUG")
    elif config.option.lahuta_info:
        logging.set_global_verbosity(LogLevel.INFO)
        print("[pytest] Lahuta logging: INFO")
    elif config.option.lahuta_error:
        logging.set_global_verbosity(LogLevel.ERROR)
        print("[pytest] Lahuta logging: ERROR")


def _resolve_data_file(filename: str) -> Path:
    here = Path(__file__).resolve()
    p = here
    for _ in range(12):
        # Try both data/ and core/data/ paths
        for base in ["core/data", "data"]:
            candidate = p.parent / base / filename if p.name == "tests" else p / base / filename
            if candidate.exists():
                return candidate
            direct = p / base / filename
            if direct.exists():
                return direct
        p = p.parent

    # Try current working directory
    for base in ["core/data", "data"]:
        alt = Path.cwd() / base / filename
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
    from lahuta import LahutaSystem

    return LahutaSystem(str(ubi_cif))


@pytest.fixture(scope="session")
def run_child() -> Callable[[str], object]:
    # Execute a Python code snippet in a subprocess and capture output
    import os
    import subprocess
    import sys
    import tempfile
    from pathlib import Path as _Path

    def _get_sanitizer_lib_path() -> str:
        """Get the TSan library path for the current platform.

        On macOS, uses clang to locate libclang_rt.tsan_osx_dynamic.dylib.
        Returns empty string if not found or on unsupported platforms.
        """
        if sys.platform != "darwin":
            return ""

        try:
            result = subprocess.run(
                ["clang", "--print-file-name=libclang_rt.tsan_osx_dynamic.dylib"],
                capture_output=True,
                text=True,
                timeout=5,
            )
            if result.returncode == 0:
                lib_path = result.stdout.strip()
                if lib_path and "/" in lib_path:
                    return lib_path
        except Exception:
            pass

        return ""

    def _inner(code: str) -> subprocess.CompletedProcess:
        env = os.environ.copy()

        sanitizer_vars = ["DYLD_INSERT_LIBRARIES", "TSAN_OPTIONS", "LSAN_OPTIONS", "UBSAN_OPTIONS", "ASAN_OPTIONS"]
        for var in sanitizer_vars:
            if var in os.environ:
                env[var] = os.environ[var]

        if "TSAN_OPTIONS" in os.environ and "DYLD_INSERT_LIBRARIES" not in env:
            lib_path = _get_sanitizer_lib_path()
            if lib_path:
                env["DYLD_INSERT_LIBRARIES"] = lib_path

        with tempfile.NamedTemporaryFile(mode="w", suffix=".py", delete=False) as f:
            f.write(code)
            temp_script = f.name

        try:
            return subprocess.run(
                [sys.executable, temp_script],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=env,
                cwd=str(_Path.cwd()),
            )
        finally:
            os.unlink(temp_script)

    return _inner


DEFAULT_LMDB = Path(__file__).resolve().parents[1] / "db_1773"


@pytest.fixture(scope="session")
def db_path() -> Path:
    env = os.getenv("LAHUTA_TEST_DB")
    if env:
        path = Path(env)
        if not path.exists():
            pytest.skip(f"LAHUTA_TEST_DB path does not exist: {path}")
        return path
    if DEFAULT_LMDB.exists():
        return DEFAULT_LMDB
    pytest.skip("Missing LMDB database (set LAHUTA_TEST_DB or ensure hum_db exists)")
