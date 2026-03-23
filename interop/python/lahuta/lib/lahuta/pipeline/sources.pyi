"""
Pipeline sources
"""
from __future__ import annotations
import typing
__all__: list[str] = ['DatabaseHandleSource', 'DatabaseSource', 'DirectorySource', 'FileListSource', 'FileSource', 'MdTrajectoriesSource', 'NmrSource', 'Source']
class DatabaseHandleSource(Source):
    def __init__(self, db: ..., database: str = '', batch: int = 1024) -> None:
        ...
class DatabaseSource(Source):
    def __init__(self, path: str, database: str = '', batch: int = 1024) -> None:
        ...
class DirectorySource(Source):
    @typing.overload
    def __init__(self, path: str, ext: str = '', recursive: bool = True, batch: int = 200) -> None:
        ...
    @typing.overload
    def __init__(self, path: str, extensions: list[str], recursive: bool = True, batch: int = 200) -> None:
        ...
class FileListSource(Source):
    def __init__(self, path: str) -> None:
        ...
class FileSource(Source):
    def __init__(self, files: list[str]) -> None:
        ...
class MdTrajectoriesSource(Source):
    def __init__(self, trajectories: typing.Any) -> None:
        """
        Create a source that streams MD trajectories (structure + XTC).
        """
class NmrSource(Source):
    def __init__(self, files: list[str]) -> None:
        """
        Create a source that streams multi-model NMR structures.
        """
class Source:
    pass
