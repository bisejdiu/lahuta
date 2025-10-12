from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Iterable, Tuple

import pytest

from lahuta import ContactProvider, InteractionType
from lahuta.pipeline import ContactTask, InMemoryPolicy, Pipeline, PipelineContext
from lahuta.sources import DirectorySource

# fmt: off
REQUIRED_ENV = "LAHUTA_CONTACTS_DIR"


def _hash_path(ctx: PipelineContext) -> str:
    return hashlib.sha256(ctx.path.encode("utf-8")).hexdigest()


@pytest.mark.skipif(
    REQUIRED_ENV not in os.environ,
    reason=f"set {REQUIRED_ENV} to run local contact parity test",
)
def test_contact_pipeline_threads_vs_processes() -> None:
    """Ensure contact outputs match between threaded and process backends for a large dataset."""
    data_dir = Path(os.environ[REQUIRED_ENV]).expanduser()
    if not data_dir.exists():
        pytest.skip(f"dataset directory does not exist: {data_dir}")

    source = DirectorySource(str(data_dir), recursive=False, extensions=[".cif.gz"], batch=64)

    pipeline = Pipeline(source)
    pipeline.params("system").is_model = True

    contact_task = ContactTask(provider=ContactProvider.GetContacts, interaction_type=InteractionType.All)
    pipeline.add_task(name="contacts",  task=contact_task, in_memory_policy=InMemoryPolicy.Keep)
    pipeline.add_task(name="path_hash", task=_hash_path, depends=["topology"], in_memory_policy=InMemoryPolicy.Keep)

    threads_result   = pipeline.run(threads=16, backend="threads")
    processes_result = pipeline.run(threads=16, backend="processes", processes=16)

    contacts_threads   = _normalize_dict_payloads(threads_result  .to_dict("contacts", columnar=False))
    contacts_processes = _normalize_dict_payloads(processes_result.to_dict("contacts", columnar=False))
    assert contacts_threads == contacts_processes

    hashes_threads   = _normalize_text_payloads(threads_result  .strings("path_hash"))
    hashes_processes = _normalize_text_payloads(processes_result.strings("path_hash"))
    assert hashes_threads == hashes_processes


def _normalize_dict_payloads(entries: Iterable[dict]) -> Tuple[str, ...]:
    dump = [json.dumps(entry, sort_keys=True, separators=(",", ":")) for entry in entries]
    return _stable_sort(dump)


def _normalize_text_payloads(entries: Iterable[str]) -> Tuple[str, ...]:
    return _stable_sort([str(entry) for entry in entries])


def _stable_sort(strings: Iterable[str]) -> Tuple[str, ...]:
    return tuple(
        sorted(
            strings,
            key=lambda s: (hashlib.sha256(s.encode("utf-8")).digest(), len(s)),
        )
    )
