from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # from lahuta.lib.lahuta.pipeline import StageManager
    from lahuta.lib import lahuta as _lib

    StageManager = _lib.pipeline.StageManager

from lahuta.lib.lahuta import TopologyComputers as TopologyComputation


class SystemParams:
    """Parameter proxy for system built-in."""

    def __init__(self, mgr: StageManager):
        self._mgr = mgr

    @property
    def is_model(self) -> bool:
        """Whether to build system from model data instead of file path."""
        return self._mgr.get_system_params()["is_model"]

    @is_model.setter
    def is_model(self, value: bool) -> None:
        self._mgr.set_system_params({"is_model": bool(value)})


class TopologyParams:
    """Parameter proxy for topology built-in."""

    def __init__(self, mgr: StageManager):
        self._mgr = mgr
        self._enabled: bool = True  # wrapper-level guard that keeps StageManager pure

    @property
    def flags(self) -> TopologyComputation:
        """Topology computation flags."""
        flags_value = self._mgr.get_topology_params()["flags"]
        # Handle uninitialized or invalid values
        if flags_value < 0:
            return TopologyComputation.All  # Default to All if uninitialized
        return TopologyComputation(flags_value)

    @flags.setter
    def flags(self, value: TopologyComputation) -> None:
        self._mgr.set_topology_params({"flags": value})

    @property
    def enabled(self) -> bool:
        """Whether the topology built-in is allowed to run.

        When False, the Python wrapper prevents adding tasks that depend on
        the "topology" builtin.
        """
        return self._enabled

    @enabled.setter
    def enabled(self, value: bool) -> None:
        self._enabled = bool(value)


__all__ = ["SystemParams", "TopologyParams"]
