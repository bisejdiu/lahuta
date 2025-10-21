from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from lahuta.lib.lahuta import ContactProvider, InteractionType, InteractionTypeSet

from .types import OutputFormat


@dataclass
class ContactTask:
    """Specification for a contacts computation."""

    provider: ContactProvider = ContactProvider.MolStar
    interaction_type: InteractionType | InteractionTypeSet | Iterable[InteractionType] = InteractionType.All
    fmt: OutputFormat = OutputFormat.BINARY

    def __post_init__(self) -> None:
        self.interaction_type = _normalize_interaction_types(self.interaction_type)


def _normalize_interaction_types(
    value: InteractionType | InteractionTypeSet | Iterable[InteractionType],
) -> InteractionTypeSet:
    if isinstance(value, InteractionTypeSet):
        return value
    if isinstance(value, InteractionType):
        return InteractionTypeSet(value)

    result = InteractionTypeSet()
    for item in value:
        if isinstance(item, InteractionTypeSet):
            result |= item
        elif isinstance(item, InteractionType):
            result |= item
        else:
            raise TypeError(
                "interaction_type must be an InteractionType, InteractionTypeSet, or an iterable of InteractionType values",
            )

    if result.empty():
        raise ValueError("interaction_type iterable must not be empty")
    return result


__all__ = ["ContactTask"]
