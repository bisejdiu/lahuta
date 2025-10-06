from __future__ import annotations

from lahuta.pipeline import PipelineContext
from lahuta.sources import FilesSource


def test_store_false_skips_context_but_still_emits(data_path) -> None:
    """Focused check: store=False does not write to TaskContext but still emits.

    Task A emits "A" but with store=False its payload is not stored under key
    'a'. Task B attempts to read ctx.get_text('a') and skips emission when None.
    """
    from lahuta.pipeline import InMemoryPolicy, Pipeline

    data = str(data_path("ubi.cif"))

    # Task A: returns a small string. We disable storing to TaskContext but keep its channel in memory
    def task_a(_: PipelineContext):
        return "A"

    # Task B: tries to read A's payload from context (should fail because store=False)
    def task_b(ctx: PipelineContext):
        # Try to read 'a' payload from context. If missing, return None to skip emission
        s = ctx.get_text("a")
        if s is None:
            return None
        return f"B:{s}"

    p = Pipeline(FilesSource([data]))

    # Add A: no storing, collect emissions in memory on channel 'a'
    p.add_task(name="a", task=task_a, depends=["system"], in_memory_policy=InMemoryPolicy.Keep, store=False)

    # Add B: depends on A and tries to read payload from A via context
    p.add_task(name="b", task=task_b, depends=["a"], in_memory_policy=InMemoryPolicy.Keep)

    out = p.run(threads=1)

    # A should have emitted exactly one value "A"
    assert "a" in out
    assert out["a"] == ["A"]

    # B should have no outputs because it could not read A's payload from context
    # (PyCallableTask returns ok=false in that case, and no emission is produced)
    assert out.get("b", []) == []


def test_store_true_allows_context_read(data_path) -> None:
    """Focused check: store=True writes payload to TaskContext for downstream reads.

    Task A emits "A" and stores it under key 'a'. Task B reads ctx.get_text('a')
    and emits a derived value.
    """
    from lahuta.pipeline import InMemoryPolicy, Pipeline

    data = str(data_path("ubi.cif"))

    def task_a(_: PipelineContext):
        return "A"

    def task_b(ctx: PipelineContext):
        s = ctx.get_text("a")
        if s is None:
            return None
        return f"B:{s}"

    p = Pipeline(FilesSource([data]))

    # Store=True: payload for 'a' is stashed into TaskContext under key "a"
    p.add_task(name="a", task=task_a, depends=["system"], in_memory_policy=InMemoryPolicy.Keep, store=True)

    # B explicitly reads from context via ctx.get_text('a')
    p.add_task(name="b", task=task_b, depends=["a"], in_memory_policy=InMemoryPolicy.Keep)

    out = p.run(threads=1)

    assert out["a"] == ["A"]
    # A returned "A" as text -> B reads text and prefixes with 'B:'
    assert out.get("b") == ["B:A"]
