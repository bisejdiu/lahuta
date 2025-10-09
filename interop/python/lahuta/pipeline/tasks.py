from __future__ import annotations

from dataclasses import dataclass

from lahuta.lib.lahuta import ContactProvider, InteractionType

from .types import OutputFormat


@dataclass
class ContactTask:
    """Specification for a contacts computation."""

    provider: ContactProvider = ContactProvider.MolStar
    interaction_type: InteractionType = InteractionType.All
    fmt: OutputFormat = OutputFormat.BINARY


__all__ = ["ContactTask"]
