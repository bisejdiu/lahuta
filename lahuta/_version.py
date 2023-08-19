"""The `version` module stores the version of lahuta."""

__all__ = ["VERSION"]

VERSION = "0.7.0"
"""The version of Lahuta."""


def get_version() -> str:
    """Get the version of Lahuta.

    Returns:
        str: The version of Lahuta.
    """
    return VERSION


def version_info() -> tuple[int, int, int]:
    """Get the version of Lahuta as a tuple of integers.

    Returns:
        tuple[int, int, int]: The version of Lahuta as a tuple of integers.
    """
    major, minor, patch = map(int, VERSION.split("."))
    return major, minor, patch
