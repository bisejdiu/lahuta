"""
Lahuta Pipeline examples:
- Build from directory, explicit files, or a filelist
- System/topology builders and dependency wiring
- Built-in contacts task with different providers and formats
- Python tasks with various inputs (path/system/topology/system_topology) and outputs (text/json)
- Channels and multiple sinks (memory, single-file, sharded files)
- Thread-safety hints and multi-threaded execution
"""

from __future__ import annotations

import os
from pathlib import Path

from lahuta import ContactProvider, InteractionType
from lahuta.pipeline import (
    ContactTask,
    FileOutput,
    InMemoryPolicy,
    OutputFormat,
    Pipeline,
    PipelineContext,
    ShardedOutput,
)

DATA_DIR = Path(__file__).resolve().parents[3] / "core" / "data"


# fmt: off
def _pick_some_files(n: int = 3) -> list[str]:
    files: list[str] = []
    for pat in ("*.cif", "*.cif.gz"):
        files.extend(str(p) for p in DATA_DIR.glob(pat))
    return files[: max(1, n)]


def ex_contacts_memory_dir() -> None:
    """Contacts over a directory, parsed JSON kept in memory"""
    p = Pipeline.from_directory(DATA_DIR, ext=".cif", recursive=False, batch=64)
    p.add_task(
        name="contacts",
        task=ContactTask(provider=ContactProvider.MolStar, interaction_type=InteractionType.All),
        in_memory_policy=InMemoryPolicy.Keep
    )
    out = p.run(threads=4)
    items = out.get("contacts", [])
    print(f"contacts (memory, parsed): {len(items)} items; first keys={list(items[0].keys()) if items and isinstance(items[0], dict) else 'n/a'}")


def ex_custom_system_stage() -> None:
    """Customize system stage: rename and depend on it"""
    p = Pipeline.from_files(_pick_some_files(1))

    # defaults to MolStar + All
    p.add_task(name="contacts", task=ContactTask(), in_memory_policy=InMemoryPolicy.Keep)
    out = p.run(threads=4)
    print(f"contacts (depends on 'init'): {len(out.get('contacts', []))}")


def ex_contacts_providers() -> None:
    """Contacts with different providers and interaction filters"""
    files = _pick_some_files(1)
    p = Pipeline.from_files(files)
    p.add_task(name="molstar_hbond", task=ContactTask(provider=ContactProvider.MolStar,  interaction_type=InteractionType.HydrogenBond), in_memory_policy=InMemoryPolicy.Keep)
    p.add_task(name="arpeggio_vdw",  task=ContactTask(provider=ContactProvider.Arpeggio, interaction_type=InteractionType.VanDerWaals),  in_memory_policy=InMemoryPolicy.Keep)
    out = p.run(threads=4)
    print("providers:",
        "molstar_hbond=", len(out.get("molstar_hbond", [])),
        "arpeggio_vdw=",  len(out.get("arpeggio_vdw",  [])),
    )


def ex_contacts_to_file_and_memory() -> None:
    """Write contacts to a single NDJSON file (and keep in memory)."""
    out_path = Path("contacts.ndjson")
    if out_path.exists():
        out_path.unlink()
    p = Pipeline.from_files(_pick_some_files(2))
    p.add_task(name="contacts",
               task=ContactTask(),
               out=[FileOutput(out_path, fmt=OutputFormat.JSON)], in_memory_policy=InMemoryPolicy.Keep)
    out = p.run(threads=4)
    files = p.file_outputs()
    print(f"contacts: memory={len(out.get('contacts', []))}, file={files.get('contacts')}")


def ex_contacts_sharded_files() -> None:
    """Sharded NDJSON outputs (no memory)."""
    out_dir = Path("contacts_shards")
    out_dir.mkdir(parents=True, exist_ok=True)

    # previous run cleanup
    for pth in out_dir.glob("part-*.ndjson"):
        pth.unlink()

    p = Pipeline.from_directory(DATA_DIR, ext=".cif")
    p.add_task(
        name="contacts",
        task=ContactTask(),
        out=[ShardedOutput(out_dir, shard_size=2, fmt=OutputFormat.JSON)],
        in_memory_policy=InMemoryPolicy.Drop
    )
    p.run(threads=4)
    shards = p.file_outputs().get("contacts", [])
    print(f"contacts shards: {len(shards)} files -> e.g., {shards[:2]}")


def ex_python_path_text() -> None:
    """Python task: path -> text"""
    def basename(ctx: PipelineContext) -> str:
        return os.path.basename(ctx.path)

    p = Pipeline.from_files(_pick_some_files(2))
    p.add_task(name="names", task=basename, in_memory_policy=InMemoryPolicy.Keep)
    out = p.run(threads=4) # can use multiple threads even with Python tasks
    print(f"python(path->text): {out.get('names', [])}")


def ex_python_path_json_parsed() -> None:
    """Python task: path -> JSON"""
    def describe(ctx: PipelineContext) -> dict:
        return {"file": os.path.basename(ctx.path), "ext": os.path.splitext(ctx.path)[1]}

    p = Pipeline.from_files(_pick_some_files(2))
    p.add_task(name="desc", task=describe, in_memory_policy=InMemoryPolicy.Keep)
    out = p.run(threads=4)
    print(f"python(path->json, parsed): {out.get('desc', [])}")


def ex_python_system_topology_input() -> None:
    """Python task: (system + topology) -> JSON"""
    def summarize(ctx: PipelineContext) -> dict:
        sys = ctx.get_system()
        top = ctx.get_topology()
        try:
            ids = list(top.get_atom_ids()) if top is not None else []
            n_atoms = int(sys.n_atoms) if sys is not None else 0
            return {"ok": True, "num_atoms": n_atoms, "first_id": int(ids[0]) if ids else None}
        except Exception:
            return {"ok": False}

    p = Pipeline.from_files(_pick_some_files(1))
    # We specify explicitly that we depend on `topology` (which itself depends on `system`)
    # This is the default behavior, and we could omit `depends` here.
    p.add_task(name="summary", task=summarize, depends=["topology"], in_memory_policy=InMemoryPolicy.Keep)
    out = p.run(threads=2)
    print(f"python(sys+topo->json): {out.get('summary', [])}")


def ex_custom_channel_multi_sinks() -> None:
    """Route a task to a custom channel and attach multiple sinks."""
    p = Pipeline.from_directory(DATA_DIR, ext=".cif", recursive=False)

    def meta(ctx: PipelineContext) -> dict[str, int | str]:
        return {"base": os.path.basename(ctx.path), "size": os.path.getsize(ctx.path)}

    p.add_task(name="meta_task", task=meta, channel="meta", in_memory_policy=InMemoryPolicy.Keep, out=[FileOutput("meta.ndjson", fmt=OutputFormat.JSON)])
    out = p.run(threads=8)
    files = p.file_outputs()
    print(f"channel 'meta': memory={len(out.get('meta', []))}, file={files.get('meta')}")


def ex_thread_safety_single_thread() -> None:
    """Thread-safety hint: single-threaded Python stage (default)"""
    def f(ctx: PipelineContext) -> str:
        return os.path.basename(ctx.path)

    p = Pipeline.from_files(_pick_some_files(5))

    # can still run with multiple threads, but this stage runs serially
    p.add_task(name="names", task=f, in_memory_policy=InMemoryPolicy.Keep, thread_safe=False)
    out = p.run(threads=8)
    print(f"thread_safe=False: processed={len(out.get('names', []))} (ran serially)")


def ex_thread_safety_multi_thread() -> None:
    """Thread-safety hint: concurrent Python stage"""
    def f(ctx: PipelineContext) -> str:
        # Pure function, safe to run concurrently
        return os.path.basename(ctx.path)

    p = Pipeline.from_files(_pick_some_files(5))
    p.add_task(name="names", task=f, in_memory_policy=InMemoryPolicy.Keep, thread_safe=True) # NOTE: thread_safe=True
    out = p.run(threads=8)
    print(f"thread_safe=True: processed={len(out.get('names', []))}")


def ex_mixed_dag() -> None:
    """Mix multiple tasks in a small DAG"""
    files = _pick_some_files(2)
    p = Pipeline.from_files(files)
    p.add_task(name="contacts", task=ContactTask(provider=ContactProvider.MolStar, interaction_type=InteractionType.All), in_memory_policy=InMemoryPolicy.Keep)

    def base(ctx: PipelineContext) -> str:
        return os.path.basename(ctx.path)

    # Independent Python stage depending only on system
    # `base` does not depend on `system` and thus it would not be included, but we explicitly add it as a dependency.
    p.add_task(name="basename", task=base, depends=["system"], in_memory_policy=InMemoryPolicy.Keep)

    def annotate(ctx: PipelineContext) -> dict[str, str | bool]:
        return {"file": os.path.basename(ctx.path), "has_contacts": True}

    # Depends on both system and contacts. Contacts depend on both system and topology, which will be included automatically.
    # This task will be called only after system, topology, and contacts are ready.
    p.add_task(name="annot", task=annotate, depends=["system", "contacts"], in_memory_policy=InMemoryPolicy.Keep)

    out = p.run(threads=4)
    print({k: len(v) for k, v in out.items()})


def ex_convenience_sinks() -> None:
    """Attach sinks via convenience methods"""
    p = Pipeline.from_files(_pick_some_files(2))
    p.add_task(name="contacts", task=ContactTask(provider=ContactProvider.MolStar, interaction_type=InteractionType.All))
    (
     p
     .to_memory       ("contacts")
     .to_files        ("contacts", path   = "contacts2.ndjson", fmt=OutputFormat.JSON)
     .to_sharded_files("contacts", out_dir= "contacts2_shards", fmt=OutputFormat.JSON, shard_size=2)
    )
    p.run(threads=4)
    files = p.file_outputs().get("contacts", [])
    print(f"convenience sinks wrote: {files}")


def ex_multiple_channels_mixed_sinks() -> None:
    """Multiple channels from multiple tasks, mixed sinks"""
    files = _pick_some_files(2)
    p = Pipeline.from_files(files)
    p.add_task(name="c_json",
               task=ContactTask(),
               out=[FileOutput("c.jsonl", fmt=OutputFormat.JSON)],
               in_memory_policy=InMemoryPolicy.Keep)
    p.add_task(name="c_sharded",
               task=ContactTask(),
               out=[ShardedOutput("c_shards", shard_size=2, fmt=OutputFormat.JSON)],
               in_memory_policy=InMemoryPolicy.Drop)

    def s(ctx: PipelineContext) -> str:
        return os.path.basename(ctx.path)

    p.add_task(name="names", task=s, in_memory_policy=InMemoryPolicy.Keep)
    out = p.run(threads=4)
    print({k: (len(v) if isinstance(v, list) else v) for k, v in out.items()})
    print("files:", p.file_outputs())


def ex_channel_fanin_contacts() -> None:
    """
    Demonstrates channel fan-in (many -> one). Two distinct contacts tasks emit to
    the same channel name "contacts", so memory results and sinks aggregate there.
    """
    files = _pick_some_files(2)
    p = Pipeline.from_files(files)

    # Two different contacts tasks, same output channel
    p.add_task(
        name="c_molstar",
        task=ContactTask(),
        channel="contacts",
        in_memory_policy=InMemoryPolicy.Keep
    )
    p.add_task(
        name="c_arpeggio",
        task=ContactTask(provider=ContactProvider.Arpeggio),
        channel="contacts",
        in_memory_policy=InMemoryPolicy.Keep
    )

    # A single sink can subscribe to the unified channel
    p.to_files("contacts", path="contacts_all.ndjson", fmt=OutputFormat.JSON)

    out = p.run(threads=4)
    files = p.file_outputs()
    print("fan-in contacts count:", len(out.get("contacts", [])))
    print("single-file sink on 'contacts':", files.get("contacts"))


def ex_channel_rename_python() -> None:
    """
    Shows how to rename a task's output stream independently of its task name.
    """

    def summarize(ctx: PipelineContext) -> dict:
        top = ctx.get_topology()
        n = len(list(top.get_atom_ids())) if top is not None else 0
        return {"atoms": n}

    p = Pipeline.from_files(_pick_some_files(1))
    p.add_task(name="py_summarize_task", task=summarize, channel="summary", in_memory_policy=InMemoryPolicy.Keep)
    out = p.run(threads=1)
    print("channels:", list(out.keys()))  # ['summary'] (not 'py_summarize_task')
    print("summary:", out.get("summary"))


def ex_channel_fanin_python_mix() -> None:
    """
    Two independent Python tasks emit into a shared channel "meta".
    Memory results and file sinks aggregate both payloads.
    """

    def base(ctx: PipelineContext) -> dict[str, str]:
        return {"file": os.path.basename(ctx.path)}

    def extra(ctx: PipelineContext) -> dict[str, int]:
        size = os.path.getsize(ctx.path)
        return {"size": int(size)}

    p = Pipeline.from_files(_pick_some_files(2))
    p.add_task(name="path_meta", task=base,  channel="meta", in_memory_policy=InMemoryPolicy.Keep)
    p.add_task(name="size_meta", task=extra, channel="meta", in_memory_policy=InMemoryPolicy.Keep)
    p.to_files("meta", path="meta.ndjson", fmt=OutputFormat.JSON)
    out = p.run(threads=1)
    print("fan-in meta items:", len(out.get("meta", [])))
    print("fan-in meta file:", p.file_outputs().get("meta"))


if __name__ == "__main__":
    ex_contacts_memory_dir()
    ex_custom_system_stage()
    ex_contacts_providers()
    ex_contacts_to_file_and_memory()
    ex_contacts_sharded_files()
    ex_python_path_text()
    ex_python_path_json_parsed()
    ex_python_system_topology_input()
    ex_custom_channel_multi_sinks()
    ex_thread_safety_single_thread()
    ex_thread_safety_multi_thread()
    ex_mixed_dag()
    ex_convenience_sinks()
    ex_multiple_channels_mixed_sinks()
