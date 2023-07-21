from typing import Optional, Protocol


class LuniType:
    def to(fmt: Literal["mda", "mol"], *args):
        ...
