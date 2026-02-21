# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     print(eval('"besian" + "sejdiu" + "@gmail.com"'))
#
from __future__ import annotations

import re
from pathlib import Path

import pytest

import lahuta
from lahuta._version import __version__ as module_version
from lahuta._version import version_info as module_version_info
from lahuta.lib import lahuta as extension


def _read_repo_version() -> str | None:
    for parent in Path(__file__).resolve().parents:
        candidate = parent / "VERSION"
        if candidate.is_file():
            return candidate.read_text(encoding="utf-8").strip()
    return None


def test_python_version_metadata_is_consistent() -> None:
    assert lahuta.__version__ == module_version
    assert lahuta.version_info.as_tuple() == module_version_info.as_tuple()

    assert getattr(extension, "__version__", None) == module_version
    assert getattr(extension, "__version_info__", None) == module_version_info.as_tuple()


def test_version_string_matches_components() -> None:
    version = module_version
    info = module_version_info

    match = re.match(r"^(\d+)\.(\d+)\.(\d+)(.*)$", version)
    assert match, f"Version '{version}' does not match expected pattern"

    major, minor, patch, suffix = match.groups()
    assert (int(major), int(minor), int(patch)) == info.as_tuple()[:3]
    assert suffix or "" == info.suffix


@pytest.mark.skipif(_read_repo_version() is None, reason="VERSION file not available alongside package")
def test_version_matches_repository_version_file() -> None:
    repo_version = _read_repo_version()
    assert repo_version is not None

    repo_match = re.match(r"^(\d+)\.(\d+)\.(\d+)(.*)$", repo_version)
    assert repo_match, f"Repository VERSION '{repo_version}' does not match expected pattern"

    repo_major, repo_minor, repo_patch, repo_suffix = repo_match.groups()
    module_major, module_minor, module_patch = module_version_info.as_tuple()[:3]

    assert (int(repo_major), int(repo_minor), int(repo_patch)) == (
        module_major,
        module_minor,
        module_patch,
    )

    normalized_suffix = module_version_info.suffix or ""
    assert bool(normalized_suffix) == bool(repo_suffix)
