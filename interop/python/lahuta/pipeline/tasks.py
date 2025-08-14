from __future__ import annotations

from dataclasses import dataclass

from lahuta.lib.lahuta import ContactProvider, InteractionType


@dataclass
class ContactTask:
    """Specification for a contacts computation."""

    provider: ContactProvider = ContactProvider.MolStar
    interaction_type: InteractionType = InteractionType.All


__all__ = ["ContactTask"]
