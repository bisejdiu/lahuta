from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Iterator, KeysView, Mapping

import orjson

from lahuta.lib import lahuta as _lib

from .types import OutputFormat


@dataclass
class _ChannelState:
    format: OutputFormat
    decoder: str | None
    data: tuple[Any, ...]

    def __str__(self) -> str:
        return f"_ChannelState(format={self.format}, decoder={self.decoder}, items={len(self.data)})"

    def __repr__(self) -> str:
        return f"_ChannelState(format={self.format}, decoder={self.decoder}, items={len(self.data)})"


class PipelineResult:
    """Encapsulates in-memory pipeline outputs with lazy decoding utilities."""

    def __init__(self, channels: Mapping[str, _ChannelState]) -> None:
        self._channels: dict[str, _ChannelState] = dict(channels)

    def __contains__(self, channel: object) -> bool:
        return channel in self._channels

    def __iter__(self) -> Iterator[str]:
        return iter(self._channels)

    def __len__(self) -> int:
        return len(self._channels)

    def __getitem__(self, channel: str) -> Any:
        state = self._require(channel)
        if state.format == OutputFormat.BINARY:
            return self.to_numpy(channel)
        if state.format == OutputFormat.JSON:
            return list(self.json(channel))
        if state.format == OutputFormat.TEXT:
            return list(self.strings(channel))
        return list(state.data)

    def get(self, channel: str, default: Any = None) -> Any:
        try:
            result = self[channel]
            # For single-item channels, return the item directly for convenience
            if isinstance(result, list) and len(result) == 1:
                return result[0]
            return result
        except KeyError:
            return default

    def channels(self) -> tuple[str, ...]:
        return tuple(self._channels)

    def keys(self) -> KeysView[str]:
        return self._channels.keys()

    def items(self) -> Iterator[tuple[str, Any]]:
        for name in self._channels:
            yield name, self[name]

    def values(self) -> Iterator[Any]:
        for name in self._channels:
            yield self[name]

    def raw(self, channel: str) -> tuple[Any, ...]:
        state = self._require(channel)
        return state.data

    def strings(self, channel: str) -> tuple[str, ...]:
        state = self._require(channel)
        if state.format == OutputFormat.BINARY:
            raise ValueError(f"Channel '{channel}' stores binary payloads. Use raw_bytes() or decode().")
        return tuple(str(item) for item in state.data)

    def json(self, channel: str, *, loads: Callable[[str], Any] | None = None) -> tuple[Any, ...]:
        state = self._require(channel)
        if state.format != OutputFormat.JSON:
            raise ValueError(f"Channel '{channel}' is not JSON formatted (has {state.format.value}).")

        if loads is not None:
            return tuple(loads(item) for item in state.data)

        payloads = tuple(str(item) for item in state.data)
        if not payloads:
            return ()

        # orjson supports batch parsing - join all items into a JSON array
        joined = "[" + ",".join(payloads) + "]"
        try:
            result = orjson.loads(joined)
        except Exception:
            # fallback to per-item parsing
            return tuple(orjson.loads(item) for item in payloads)

        if isinstance(result, list):
            return tuple(result)
        elif isinstance(result, tuple):
            return result
        else:
            return (result,)

    def _contacts_numpy(self, channel: str) -> tuple[Any, ...]:
        state = self._require(channel)
        if state.decoder != "contacts":
            raise ValueError(f"Channel '{channel}' is not registered for contacts decoding.")
        if state.format != OutputFormat.BINARY:
            raise ValueError(f"Channel '{channel}' stores {state.format.value} data, not binary contacts payloads.")

        return tuple(_lib.decode_contacts_binary(payload) for payload in state.data)

    def to_numpy(self, channel: str | None = None) -> tuple[Any, ...]:
        """Return columnar (NumPy-backed) contact payloads for the selected channel."""
        name = self._resolve_channel(channel)
        state = self._require(name)
        if state.decoder != "contacts" or state.format != OutputFormat.BINARY:
            raise ValueError(
                f"Channel '{name}' does not expose a NumPy decoder. "
                "Only contact channels stored in binary format support to_numpy()."
            )
        return self._contacts_numpy(name)

    def to_dict(self, channel: str | None = None, columnar: bool = False) -> list[Any]:
        """
        Return pure-Python contact or JSON/text payloads for the selected channel.

        Args:
            channel: Channel name (optional if only one channel)
            columnar: If True, returns dict-of-lists format for better performance.
        """
        name = self._resolve_channel(channel)
        state = self._require(name)

        if state.format == OutputFormat.JSON:
            return list(self.json(name))
        if state.format == OutputFormat.TEXT:
            return [str(item) for item in state.data]
        if state.format == OutputFormat.BINARY:
            if state.decoder != "contacts":
                raise ValueError(
                    f"Channel '{name}' stores binary data without a registered decoder. "
                    "Use raw_bytes() or provide a custom decoder."
                )

            # Use parallel batch decoding for large datasets (>= 50 items)
            # For smaller datasets, sequential is faster due to thread overhead
            if len(state.data) >= 50 and hasattr(_lib, "decode_contacts_batch_parallel"):
                # Parallel batch processing - releases GIL and uses multiple cores
                return _lib.decode_contacts_batch_parallel(list(state.data), columnar)
            else:
                # Sequential processing for small batches
                if columnar:
                    return [_lib.decode_contacts_binary_columnar(payload) for payload in state.data]
                else:
                    return [_lib.decode_contacts_binary_direct(payload) for payload in state.data]

        raise ValueError(f"Unsupported channel format {state.format} for to_dict().")

    def _require(self, channel: str) -> _ChannelState:
        try:
            return self._channels[channel]
        except KeyError as exc:
            raise KeyError(f"Channel '{channel}' not in result. Available: {tuple(self._channels)}") from exc

    def _resolve_channel(self, channel: str | None) -> str:
        if channel is not None:
            return channel
        if len(self._channels) == 1:
            return next(iter(self._channels))
        raise ValueError("Multiple channels present. Please specify which channel to use (e.g. .to_dict('contacts')).")


__all__ = ["PipelineResult"]
