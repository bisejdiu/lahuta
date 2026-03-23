"""
LMDB integration
"""
from __future__ import annotations
import lahuta.lib.lahuta
__all__: list[str] = ['Database']
class Database:
    def __init__(self, path: str, max_size_gb: int = 500) -> None:
        ...
    def get_model(self, key: str) -> lahuta.lib.lahuta.LahutaSystem:
        """
        Create and return a LahutaSystem built from the stored model
        """
    def get_raw(self, key: str) -> bytes:
        ...
    def keys(self) -> list[str]:
        """
        Return all keys in the database (eager list)
        """
