from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from lahuta.lib.lahuta import ContactProvider, InteractionType, InteractionTypeSet

from .._interaction_utils import normalize_interaction_selection
from .types import OutputFormat


@dataclass
class ContactTask:
    """Specification for a contacts computation."""

    provider: ContactProvider = ContactProvider.MolStar
    interaction_type: InteractionType | InteractionTypeSet | Iterable[InteractionType] = InteractionType.All
    fmt: OutputFormat = OutputFormat.BINARY

    def __post_init__(self) -> None:
        self.interaction_type = normalize_interaction_selection(self.interaction_type, arg_name="interaction_type")


__all__ = ["ContactTask"]
