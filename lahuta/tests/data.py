from .base import BaseFile

__all__ = ["X2", "Rhodopsin", "DNABound"]


class X2(BaseFile):
    FILE_NAME = "1KX2"

    def __init__(self) -> None:
        super().__init__(pdb=True)


class Rhodopsin(BaseFile):
    FILE_NAME = "1GZM"


class DNABound(BaseFile):
    FILE_NAME = "3Q2Y"
