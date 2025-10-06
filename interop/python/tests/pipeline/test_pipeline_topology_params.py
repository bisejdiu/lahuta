from __future__ import annotations

import tempfile
from pathlib import Path
from typing import TypedDict

import pytest

from lahuta.lib.lahuta import TopologyComputers as TopologyComputation
from lahuta.pipeline import Pipeline
from lahuta.sources import DirectorySource, FileListSource, FileSource


class ProbeRec(TypedDict):
    has_system: bool
    has_topology: bool


class InspectRec(TypedDict):
    neighbors: bool
    bonds: bool


class InspectRecFull(InspectRec):
    atom_typing: bool


class NeedsRec(TypedDict):
    has_system: bool


DATA_FILE = Path(__file__).resolve().parents[4] / "core" / "data" / "1kx2_small.cif"


def create_test_files_in_temp_dir(temp_dir: Path, count: int = 2) -> list[str]:
    """Create temporary test files for testing different source types."""

    test_files = []
    for i in range(count):
        test_file = temp_dir / f"test_{i}.cif"
        # Copy the source file
        with open(DATA_FILE, "r") as src, open(test_file, "w") as dst:
            dst.write(src.read())
        test_files.append(str(test_file))

    return test_files


def create_file_list(temp_dir: Path, files: list[str]) -> str:
    """Create a file list for FileListSource testing."""
    file_list = temp_dir / "file_list.txt"
    with open(file_list, "w") as f:
        for file_path in files:
            f.write(f"{file_path}\n")
    return str(file_list)


def test_builtins_not_run_when_not_required(single_test_file: str) -> None:
    """If no task depends on built-ins, they are not executed."""

    p = Pipeline(FileSource([single_test_file]))

    def probe(ctx) -> ProbeRec:
        sys = ctx.get_system()
        top = ctx.get_topology()
        return {"has_system": bool(sys), "has_topology": bool(top)}

    # Explicitly override auto-discovery, which means no built-ins will be injected
    p.add_task(name="probe", task=probe, depends=[])

    class Schema(TypedDict):
        probe: list[ProbeRec]

    out = p.run_typed(Schema, threads=1)
    assert "probe" in out and len(out["probe"]) == 1
    rec = out["probe"][0]
    assert rec["has_system"] is False
    assert rec["has_topology"] is False


def test_topology_flags_neighbors_only(minimal_test_files: list[str]) -> None:
    """Users can restrict topology via flags; only requested computations are enabled."""

    p = Pipeline(FileSource(minimal_test_files))
    p.params("topology").flags = TopologyComputation.Neighbors

    def inspect(ctx) -> InspectRec:
        top = ctx.get_topology()
        assert top is not None
        return {
            "neighbors": top.is_computation_enabled(TopologyComputation.Neighbors),
            "bonds": top.is_computation_enabled(TopologyComputation.Bonds),
        }

    p.add_task(name="inspect", task=inspect, depends=["topology"])  # requires topology

    class Schema(TypedDict):
        inspect: list[InspectRec]

    out = p.run_typed(Schema, threads=1)
    assert "inspect" in out and len(out["inspect"]) == len(minimal_test_files)
    for rec in out["inspect"]:
        assert rec["neighbors"] is True
        assert rec["bonds"] is False


def test_disable_topology_guard(minimal_test_files: list[str]) -> None:
    """Disabling topology prevents adding tasks that require it (wrapper-level guard)."""

    p = Pipeline(FileSource(minimal_test_files))
    p.params("topology").enabled = False

    # Explicit dependency on 'topology' should be rejected
    with pytest.raises(ValueError):
        p.add_task(name="inspect", task=lambda ctx: None, depends=["topology"])  # noqa: E731

    # ContactsTask requires topology; should also be rejected
    from lahuta.lib.lahuta import ContactProvider, InteractionType
    from lahuta.pipeline import ContactTask

    with pytest.raises(ValueError):
        p.add_task(name="contacts", task=ContactTask(ContactProvider.MolStar, InteractionType.All))


def test_implicit_system_dependency(single_test_file: str) -> None:
    """Depending on 'system' works without registering a task; engine provides it on-demand."""

    p = Pipeline(FileSource([single_test_file]))

    def needs_system(ctx) -> NeedsRec:
        return {"has_system": bool(ctx.get_system())}

    p.add_task(name="needs_system", task=needs_system, depends=["system"])  # implicit built-in

    class Schema(TypedDict):
        needs_system: list[NeedsRec]

    out = p.run_typed(Schema, threads=1)
    assert out["needs_system"][0]["has_system"] is True


def test_toggle_topology_flags_invalidation_vector_source(single_test_file: str) -> None:
    """Test VectorSource (from_files) - toggling flags between runs should update enabled computations."""

    p = Pipeline(FileSource([single_test_file]))

    def inspect(ctx) -> InspectRec:
        top = ctx.get_topology()
        assert top is not None
        return {
            "neighbors": top.is_computation_enabled(TopologyComputation.Neighbors),
            "bonds": top.is_computation_enabled(TopologyComputation.Bonds),
        }

    p.add_task(name="inspect", task=inspect, depends=["topology"])  # requires topology

    # First run: neighbors only
    p.params("topology").flags = TopologyComputation.Neighbors

    class Schema(TypedDict):
        inspect: list[InspectRec]

    out1 = p.run_typed(Schema, threads=1)
    rec1 = out1["inspect"][0]
    assert rec1["neighbors"] is True
    assert rec1["bonds"] is False

    # Second run: bonds (engine should reflect new flags)
    p.params("topology").flags = TopologyComputation.Bonds
    out2 = p.run_typed(Schema, threads=1)
    # Memory sinks are cleared between runs; for multi-item sources,
    # choosing the last record is a simple way to select the most recent item.
    rec2 = out2["inspect"][-1]
    assert rec2["bonds"] is True
    # We do not assert neighbors explicitly here because engine may auto-enable
    # neighbors as a dependency of bonds at execution time.


def test_toggle_topology_flags_invalidation_directory_source() -> None:
    """Test DirectorySource (from_directory) with multiple runs and parameter changes."""
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        _ = create_test_files_in_temp_dir(temp_path, count=2)

        source = DirectorySource(str(temp_path), extensions=[".cif"], recursive=False)
        p = Pipeline(source)

        def inspect(ctx) -> InspectRecFull:
            top = ctx.get_topology()
            assert top is not None
            return {
                "neighbors": top.is_computation_enabled(TopologyComputation.Neighbors),
                "bonds": top.is_computation_enabled(TopologyComputation.Bonds),
                "atom_typing": top.is_computation_enabled(TopologyComputation.AtomTyping),
            }

        p.add_task(name="inspect", task=inspect, depends=["topology"])

        # First run: Neighbors only
        p.params("topology").flags = TopologyComputation.Neighbors

        class Schema(TypedDict):
            inspect: list[InspectRecFull]

        out1 = p.run_typed(Schema, threads=1)
        assert "inspect" in out1 and out1["inspect"], "DirectorySource: No results from first run"
        # Should have 2 results (one per file)
        assert len(out1["inspect"]) == 2, f"DirectorySource: Expected 2 results, got {len(out1['inspect'])}"

        for rec in out1["inspect"]:
            assert rec["neighbors"] is True
            assert rec["bonds"] is False
            assert rec["atom_typing"] is False

        # Second run: Bonds only
        p.params("topology").flags = TopologyComputation.Bonds
        out2 = p.run_typed(Schema, threads=1)
        assert "inspect" in out2 and out2["inspect"], "DirectorySource: No results from second run (source exhausted?)"

        # Check latest results
        for rec in out2["inspect"][-2:]:  # Get last 2 results
            assert rec["bonds"] is True

        # Third run: All flags
        p.params("topology").flags = TopologyComputation.All
        out3 = p.run_typed(Schema, threads=1)
        assert "inspect" in out3 and out3["inspect"], "DirectorySource: No results from third run (source exhausted?)"

        for rec in out3["inspect"][-2:]:  # Get last 2 results
            assert rec["neighbors"] is True
            assert rec["bonds"] is True
            assert rec["atom_typing"] is True


def test_toggle_topology_flags_invalidation_filelist_source() -> None:
    """Test FileListSource (from_filelist) with multiple runs and parameter changes."""
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        test_files = create_test_files_in_temp_dir(temp_path, count=2)
        file_list_path = create_file_list(temp_path, test_files)

        p = Pipeline(FileListSource(file_list_path))

        def inspect(ctx) -> InspectRecFull:
            top = ctx.get_topology()
            assert top is not None
            return {
                "neighbors": top.is_computation_enabled(TopologyComputation.Neighbors),
                "bonds": top.is_computation_enabled(TopologyComputation.Bonds),
                "atom_typing": top.is_computation_enabled(TopologyComputation.AtomTyping),
            }

        p.add_task(name="inspect", task=inspect, depends=["topology"])

        # First run: Neighbors only
        p.params("topology").flags = TopologyComputation.Neighbors

        class Schema(TypedDict):
            inspect: list[InspectRecFull]

        out1 = p.run_typed(Schema, threads=1)
        assert "inspect" in out1 and out1["inspect"], "FileListSource: No results from first run"
        assert len(out1["inspect"]) == 2, f"FileListSource: Expected 2 results, got {len(out1['inspect'])}"

        for rec in out1["inspect"]:
            assert rec["neighbors"] is True
            assert rec["bonds"] is False
            assert rec["atom_typing"] is False

        # Second run: Bonds only
        p.params("topology").flags = TopologyComputation.Bonds
        out2 = p.run_typed(Schema, threads=1)
        assert "inspect" in out2 and out2["inspect"], "FileListSource: No results from second run (source exhausted?)"

        for rec in out2["inspect"][-2:]:  # Get last 2 results
            assert rec["bonds"] is True

        # Third run: All flags
        p.params("topology").flags = TopologyComputation.All
        out3 = p.run_typed(Schema, threads=1)
        assert "inspect" in out3 and out3["inspect"], "FileListSource: No results from third run (source exhausted?)"

        for rec in out3["inspect"][-2:]:  # Get last 2 results
            assert rec["neighbors"] is True
            assert rec["bonds"] is True
            assert rec["atom_typing"] is True


@pytest.mark.parametrize("source_type", ["vector", "directory", "filelist"])
def test_all_sources_multi_run_parameter_changes(source_type: str) -> None:
    """Parametrized test for all source types to verify multi-run parameter invalidation works."""
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        test_files = create_test_files_in_temp_dir(temp_path, count=1)  # Single file for simplicity

        # Create pipeline based on source type
        if source_type == "vector":
            p = Pipeline(FileSource(test_files))
        elif source_type == "directory":
            p = Pipeline(DirectorySource(str(temp_path), extensions=[".cif"], recursive=False))
        elif source_type == "filelist":
            file_list = create_file_list(temp_path, test_files)
            p = Pipeline(FileListSource(file_list))
        else:
            pytest.skip(f"Unknown source type: {source_type}")

        def inspect(ctx) -> InspectRec:
            top = ctx.get_topology()
            assert top is not None
            return {
                "neighbors": top.is_computation_enabled(TopologyComputation.Neighbors),
                "bonds": top.is_computation_enabled(TopologyComputation.Bonds),
            }

        p.add_task(name="inspect", task=inspect, depends=["topology"])

        # First run: neighbors only
        p.params("topology").flags = TopologyComputation.Neighbors

        class Schema(TypedDict):
            inspect: list[InspectRec]

        out1 = p.run_typed(Schema, threads=1)
        assert "inspect" in out1 and out1["inspect"], f"{source_type}: No results from first run"
        rec1 = out1["inspect"][0]
        assert rec1["neighbors"] is True
        assert rec1["bonds"] is False

        # Second run: bonds - this will fail for sources without proper reset() capability
        p.params("topology").flags = TopologyComputation.Bonds
        out2 = p.run_typed(Schema, threads=1)

        # This assertion will fail for broken sources
        assert "inspect" in out2 and out2["inspect"], (
            f"{source_type}: No results from second run - source is not reusable!"
        )
        rec2 = out2["inspect"][-1]
        assert rec2["bonds"] is True, f"{source_type}: Flags not updated on second run"
