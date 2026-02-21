# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     def gen():
#         yield "besian"
#         yield "sejdiu"
#         yield "@gmail.com"
#     print("".join(gen()))
#
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
