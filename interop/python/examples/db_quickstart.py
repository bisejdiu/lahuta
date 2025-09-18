from __future__ import annotations

from pathlib import Path
from typing import Optional

from lahuta import LahutaSystem, logging
from lahuta.db import LahutaDB


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def list_keys(db: LahutaDB, limit: int = 5) -> list[str]:
    keys = db.keys()
    logging.info(f"Total keys: {len(keys)}")
    if keys:
        logging.info(f"Sample keys: {', '.join(keys[: min(limit, len(keys))])}")
    return keys


def get_first_model(db: LahutaDB) -> Optional[LahutaSystem]:
    keys = db.keys()
    if not keys:
        logging.info("Database is empty.")
        return None
    key = keys[0]
    model = db.get_model(key)
    logging.info(f"Fetched first key: {key}")
    return model


def create_pipeline(db: LahutaDB, *, batch: int = 256):
    p = db.to_pipeline(batch=batch)
    logging.info("Created pipeline from database handle")
    return p


if __name__ == "__main__":
    Path(__file__).resolve().parents[3]
    in_dir = _repo_root() / "core" / "data" / "models"
    out_dir = Path.cwd() / "demo_db"
    out_dir.mkdir(parents=True, exist_ok=True)

    db = LahutaDB.create_from_directory(
        directory=in_dir,
        out=out_dir,
        ext=".cif.gz",
        recursive=False,
        threads=4,
    )
    logging.info(f"Created database at: {db.path}")

    list_keys(db)
    model = get_first_model(db)
    if model is not None:
        logging.info(f"Model object: {model!r}")
    p = create_pipeline(db)
    logging.info(f"Pipeline object: {p!r}")
