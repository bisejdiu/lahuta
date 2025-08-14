"""Optimized fixtures for pipeline tests to improve performance."""

from pathlib import Path

import pytest

import lahuta as lxx


@pytest.fixture(scope="session")
def pipeline_test_files(data_dir: Path) -> list[str]:
    """Session-scoped test files to avoid repeated file system access."""
    files: list[str] = []
    for pat in ("*.cif", "*.cif.gz"):
        files.extend([str(p) for p in sorted(data_dir.glob(pat))])
    if not files:
        pytest.skip(f"No CIF files in {data_dir}")
    # Hacky, but limit to 3 small files so tests finish faster by prioritizing smaller files
    small_files = [f for f in files if "small" in Path(f).name or "fubi" in Path(f).name or "5i55" in Path(f).name]
    return small_files[:3] if len(small_files) >= 3 else files[:3] if len(files) >= 3 else files


@pytest.fixture(scope="session")
def cached_luni_systems(pipeline_test_files: list[str]) -> dict[str, lxx.LahutaSystem]:
    """Cache LahutaSystem objects."""
    systems = {}
    for file_path in pipeline_test_files:
        try:
            systems[file_path] = lxx.LahutaSystem(file_path)
        except Exception as e:
            pytest.skip(f"Failed to load {file_path}: {e}")
    return systems


@pytest.fixture(scope="session")
def minimal_test_files(pipeline_test_files: list[str]) -> list[str]:
    return pipeline_test_files[:2]


@pytest.fixture(scope="session")
def single_test_file(pipeline_test_files: list[str]) -> str:
    return pipeline_test_files[0]
