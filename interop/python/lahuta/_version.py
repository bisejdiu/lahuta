from __future__ import annotations

import re
from dataclasses import dataclass
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Iterable

__all__ = ["VersionInfo", "__version__", "version_info"]


@dataclass(frozen=True)
class VersionInfo:
    """Lahuta version information."""

    major: int
    minor: int
    patch: int
    suffix: str

    def as_tuple(self) -> tuple[int, int, int, str]:
        return (self.major, self.minor, self.patch, self.suffix)


_VERSION_PATTERN = re.compile(r"^(?P<major>\d+)\.(?P<minor>\d+)\.(?P<patch>\d+)(?P<suffix>.*)$")

_CMAKE_VERSION = "@LAHUTA_VERSION_FULL@"
_CMAKE_CONFIGURED = "@LAHUTA_CMAKE_CONFIGURED@"


def _candidate_version_files() -> Iterable[Path]:
    for directory in Path(__file__).resolve().parents:
        candidate = directory / "VERSION"
        if candidate.is_file():
            yield candidate


def _read_version_from_package() -> str | None:
    try:
        return importlib_metadata.version("lahuta")
    except importlib_metadata.PackageNotFoundError:
        return None


def _read_version_from_tree() -> str | None:
    for candidate in _candidate_version_files():
        text = candidate.read_text(encoding="utf-8").strip()
        if text:
            return text
    return None


def _resolve_version() -> str:
    if _CMAKE_CONFIGURED == "TRUE" and _CMAKE_VERSION and not _CMAKE_VERSION.startswith("@"):
        return _CMAKE_VERSION

    # Try package metadata or VERSION file
    for supplier in (_read_version_from_package, _read_version_from_tree):
        value = supplier()
        if value:
            return value
    return "0.0.0"


__version__ = _resolve_version()

_match = _VERSION_PATTERN.match(__version__)
if _match:
    version_info = VersionInfo(
        major=int(_match.group("major")),
        minor=int(_match.group("minor")),
        patch=int(_match.group("patch")),
        suffix=_match.group("suffix") or "",
    )
else:
    version_info = VersionInfo(0, 0, 0, __version__)
