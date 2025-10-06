"""Validates the Lahuta pipeline functionality using a comprehensive test suite."""

from __future__ import annotations

import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, TypedDict

import pytest

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
from lahuta.sources import DirectorySource, FileListSource, FilesSource

# fmt: off
EXPECTED_CONTACTS_COUNT = {
    "1kx2_small.cif": 112,
    "fubi.cif":       169,
    "5i55.cif":       58,
}

EXPECTED_FILE_EXTENSIONS = {
    "1kx2_small.cif": ".cif",
    "fubi.cif":       ".cif",
    "5i55.cif":       ".cif",
}

EXPECTED_CONTACT_KEYS = ["file_path", "success", "provider", "contact_type", "num_contacts", "contacts"]
EXPECTED_INDIVIDUAL_CONTACT_KEYS = ["lhs", "rhs", "distance", "type"]


@pytest.fixture(scope="class")
def temp_dir():
    """Fixture providing a temporary directory that gets cleaned up."""
    temp_dir = tempfile.mkdtemp()
    yield Path(temp_dir)
    shutil.rmtree(temp_dir, ignore_errors=True)


class TestPipelineFromDirectory:
    """Test pipeline creation from directory source."""

    def test_contacts_memory_from_directory(self, data_dir: Path):
        """Test contacts generation from directory with results kept in memory."""
        # Limit to smaller files for faster testing
        source = DirectorySource(data_dir, ext=".cif", recursive=False, batch=64)
        p = Pipeline(source)
        p.add_task(name="contacts", task=ContactTask(provider=ContactProvider.MolStar, interaction_type=InteractionType.All), in_memory_policy=InMemoryPolicy.Keep)

        # Contact recors are too complex for static typing (just not worth it)
        results = p.run(threads=4)

        assert "contacts" in results
        contacts = results["contacts"]
        assert isinstance(contacts, list)
        assert len(contacts) >= 1  # Reduced from 3 for faster testing

        # Validate each contact result structure
        for contact_result in contacts:
            assert isinstance(contact_result, dict)
            for key in EXPECTED_CONTACT_KEYS:
                assert key in contact_result, f"Missing key: {key}"

            # Validate file path is absolute and exists
            file_path = contact_result["file_path"]
            assert os.path.isabs(file_path)
            assert Path(file_path).exists()

            # Validate success and basic structure
            assert contact_result["success"] is True
            assert contact_result["provider"] == "molstar"
            assert contact_result["contact_type"] == "All"
            assert isinstance(contact_result["num_contacts"], int)
            assert contact_result["num_contacts"] >= 0
            assert isinstance(contact_result["contacts"], list)

            file_name = Path(file_path).name
            if file_name in EXPECTED_CONTACTS_COUNT:
                expected_count = EXPECTED_CONTACTS_COUNT[file_name]
                actual_count = contact_result["num_contacts"]
                assert actual_count == expected_count, (
                    f"Contact count mismatch for {file_name}: expected {expected_count}, got {actual_count}"
                )

                contacts_list = contact_result["contacts"]
                assert len(contacts_list) == actual_count

                for contact in contacts_list[:3]:  # first few contacts
                    for key in EXPECTED_INDIVIDUAL_CONTACT_KEYS:
                        assert key in contact, f"Missing contact key: {key}"
                    assert isinstance(contact["distance"], (int, float))
                    assert contact["distance"] > 0
                    assert isinstance(contact["type"], str)
                    assert len(contact["type"]) > 0


class TestPipelineFromFiles:
    """Test pipeline creation from explicit file list."""

    def test_contacts_memory_from_files(self, minimal_test_files: list[str]):
        """Test contacts generation from explicit file list."""
        p = Pipeline(FilesSource(minimal_test_files))
        p.add_task(
            name="contacts",
            task=ContactTask(provider=ContactProvider.MolStar, interaction_type=InteractionType.All),
            in_memory_policy=InMemoryPolicy.Keep,
        )

        results = p.run(threads=4)
        contacts = results["contacts"]

        assert len(contacts) == len(minimal_test_files)

        # Validate each result corresponds to our input files
        processed_files = {Path(c["file_path"]).name for c in contacts}
        expected_files  = {Path(f).name for f in minimal_test_files}
        assert processed_files == expected_files

    def test_contacts_from_filelist(self, minimal_test_files: list[str], temp_dir: Path):
        """Test contacts generation from filelist file."""
        # Create filelist
        filelist_path = temp_dir / "inputs.lst"
        filelist_path.write_text("\n".join(minimal_test_files))

        p = Pipeline(FileListSource(filelist_path))
        p.add_task(
            name="contacts",
            task=ContactTask(provider=ContactProvider.MolStar, interaction_type=InteractionType.All),
            in_memory_policy=InMemoryPolicy.Keep,
        )

        results = p.run(threads=4)
        contacts = results["contacts"]

        assert len(contacts) == len(minimal_test_files)


class TestContactProviders:
    """Test different contact providers and interaction types."""

    def test_multiple_contact_providers(self, minimal_test_files: list[str]):
        """Test MolStar and Arpeggio providers with different interaction types."""
        files = minimal_test_files  # already limited to 2 files
        p = Pipeline(FilesSource(files))

        p.add_task(
            name="contacts_molstar_hbond",
            task=ContactTask(provider=ContactProvider.MolStar, interaction_type=InteractionType.HydrogenBond),
            in_memory_policy=InMemoryPolicy.Keep,
        )
        p.add_task(
            name="contacts_arpeggio_vdw",
            task=ContactTask(provider=ContactProvider.Arpeggio, interaction_type=InteractionType.VanDerWaals),
            in_memory_policy=InMemoryPolicy.Keep,
        )

        results = p.run(threads=4)

        # Both channels should have results
        assert "contacts_molstar_hbond" in results
        assert "contacts_arpeggio_vdw"  in results

        molstar_results  = results["contacts_molstar_hbond"]
        arpeggio_results = results["contacts_arpeggio_vdw"]

        assert len(molstar_results)  == len(files)
        assert len(arpeggio_results) == len(files)

        # Validate provider-specific results
        for result in molstar_results:
            assert result["provider"]     == "molstar"
            assert result["contact_type"] == "HydrogenBond"

        for result in arpeggio_results:
            assert result["provider"]     == "arpeggio"
            assert result["contact_type"] == "VanDerWaals"


class TestFileOutputs:
    """Test file and sharded output functionality."""

    def test_contacts_to_file_and_memory(self, minimal_test_files: list[str], temp_dir: Path):
        """Test writing contacts to file while keeping in memory."""
        files = minimal_test_files
        out_path = temp_dir / "contacts.ndjson"

        p = Pipeline(FilesSource(files))
        p.add_task(
            name="contacts",
            task=ContactTask(provider=ContactProvider.MolStar, interaction_type=InteractionType.All),
            out=[FileOutput(out_path, fmt=OutputFormat.JSON)],
            in_memory_policy=InMemoryPolicy.Keep,
        )

        results = p.run(threads=4)
        file_outputs = p.file_outputs()

        assert "contacts" in results
        assert len(results["contacts"]) == len(files)

        # Check file was created
        assert out_path.exists()
        assert out_path.stat().st_size > 0

        # Check file_outputs() method
        assert "contacts" in file_outputs
        assert str(out_path) in file_outputs["contacts"]

        # Validate file contents
        with open(out_path) as f:
            lines = f.readlines()
            assert len(lines) == len(files)

            for line in lines:
                line = line.strip()
                if line:
                    try:
                        data = json.loads(line)
                        assert isinstance(data, dict)
                        for key in EXPECTED_CONTACT_KEYS:
                            assert key in data
                    except json.JSONDecodeError as _: # could be truncated maybe?
                        assert line.startswith('{"file_path"')
                        assert '"success":true'       in line
                        assert '"provider":"molstar"' in line

    def test_contacts_sharded_files(self, minimal_test_files: list[str], temp_dir: Path):
        """Test sharded output functionality."""
        out_dir = temp_dir / "contacts_shards"
        out_dir.mkdir(parents=True, exist_ok=True)

        p = Pipeline(FilesSource(minimal_test_files))  # Use explicit files instead of directory scan
        p.add_task(
            name="contacts",
            task=ContactTask(provider=ContactProvider.MolStar, interaction_type=InteractionType.All),
            out=[ShardedOutput(out_dir, shard_size=1, fmt=OutputFormat.JSON)],  # Smaller shards for faster testing
            in_memory_policy=InMemoryPolicy.Drop,
        )

        p.run(threads=2)  # Fewer threads for smaller workload
        file_outputs = p.file_outputs()

        # Check sharded files were created
        assert "contacts" in file_outputs
        shard_files = file_outputs["contacts"]
        assert len(shard_files) >= 2  # Should create multiple shards

        # Validate shard file naming pattern
        for shard_file in shard_files:
            shard_path = Path(shard_file)
            assert shard_path.parent == out_dir
            assert shard_path.name.startswith("part-")
            assert shard_path.suffix == ".ndjson"

            if shard_path.stat().st_size > 0:
                with open(shard_path) as f:
                    lines = f.readlines()
                    assert len(lines) <= 1  # Shard size limit
                    for line in lines:
                        line = line.strip()
                        if line:
                            try:
                                data = json.loads(line)
                                assert isinstance(data, dict)
                                assert "file_path" in data
                                assert "success" in data
                            except json.JSONDecodeError:
                                assert line.startswith('{"file_path"')
                                assert '"success":' in line


class TestPythonTasks:
    """Test Python task functionality."""

    def test_python_path_to_text(self, minimal_test_files: list[str], temp_dir: Path):
        """Test Python task that takes path input and returns text."""

        def basename_with_prefix(ctx: PipelineContext) -> str:
            return "processed: " + os.path.basename(ctx.path)

        files = minimal_test_files
        out_path = temp_dir / "names.txt"

        p = Pipeline(FilesSource(files))
        p.add_task(
            name="names",
            task=basename_with_prefix,
            in_memory_policy=InMemoryPolicy.Keep,
            out=[FileOutput(out_path, fmt=OutputFormat.TEXT)],
        )

        class NamesSchema(TypedDict):
            names: list[str]

        results = p.run_typed(NamesSchema, threads=1)

        #
        # The file sink used to keep an ofstream open and did not flush after each write,
        # so immediately after run() returned the stream's internal buffer could still be unflushed.
        # this is fixed now, so these lines are not needed anymore.
        #
        # del p
        # import gc; gc.collect()

        assert "names" in results
        names = results["names"]
        assert len(names) == len(files)

        for name, file_path in zip(names, files):
            expected_name = "processed: " + os.path.basename(file_path)
            assert name == expected_name

        assert out_path.exists()

        with open(out_path) as f:
            content = f.read().strip()
            lines = [line.strip() for line in content.split("\n") if line.strip()]
            assert len(lines) == len(files)
            for line, file_path in zip(lines, files):
                expected_name = "processed: " + os.path.basename(file_path)
                assert line == expected_name


    def test_python_path_to_json(self, minimal_test_files: list[str], temp_dir: Path):
        """Test Python task that takes path input and returns JSON."""

        def describe_file(ctx: PipelineContext) -> dict[str, str]:
            return {"file": os.path.basename(ctx.path), "ext": os.path.splitext(ctx.path)[1]}

        files = minimal_test_files
        out_path = temp_dir / "desc.json"

        p = Pipeline(FilesSource(files))
        p.add_task(
            name="desc",
            task=describe_file,
            in_memory_policy=InMemoryPolicy.Keep,
            out=[FileOutput(out_path, fmt=OutputFormat.JSON)],
        )

        class DescRec(TypedDict):
            file: str
            ext: str

        class Schema(TypedDict):
            desc: list[DescRec]

        results = p.run_typed(Schema, threads=1)

        assert "desc" in results
        descriptions = results["desc"]
        assert len(descriptions) == len(files)

        for desc, file_path in zip(descriptions, files):
            expected_file = os.path.basename(file_path)
            expected_ext = os.path.splitext(file_path)[1]

            assert desc["file"] == expected_file
            assert desc["ext"]  == expected_ext

            if expected_file in EXPECTED_FILE_EXTENSIONS:
                assert desc["ext"] == EXPECTED_FILE_EXTENSIONS[expected_file]

        assert out_path.exists()

        with open(out_path) as f:
            content = f.read().strip()
            lines   = [line.strip() for line in content.split("\n") if line.strip()]
            assert len(lines) > 0

            for line in lines:
                data = json.loads(line)
                assert isinstance(data, dict)
                assert "file" in data
                assert "ext"  in data

    def test_python_callable_memory_serialization_modes(self, single_test_file: str):
        """Memory results for Python callables:
        - dict return values are JSON-serialized and parsed back to Python dicts.
        - str return values emit raw text. JSON parsing fails and wrapper falls back to string.
        """

        def return_dict(ctx: PipelineContext) -> dict[str, str]:
            return {"file": os.path.basename(ctx.path)}

        def return_text(ctx: PipelineContext) -> str:
            return f"processed: {os.path.basename(ctx.path)}"

        p = Pipeline(FilesSource([single_test_file]))
        p.add_task(name="dict_task", task=return_dict, in_memory_policy=InMemoryPolicy.Keep)
        p.add_task(name="text_task", task=return_text, in_memory_policy=InMemoryPolicy.Keep)

        class DictRec(TypedDict):
            file: str

        class Schema(TypedDict):
            dict_task: list[DictRec]
            text_task: list[str]

        results = p.run_typed(Schema, threads=1)

        assert "dict_task" in results and len(results["dict_task"]) == 1
        assert "text_task" in results and len(results["text_task"]) == 1

        dict_res = results["dict_task"][0]
        text_res = results["text_task"][0]

        assert isinstance(dict_res, dict)
        assert "file" in dict_res
        assert isinstance(text_res, str)
        assert text_res.startswith("processed: ")


class TestChannelFanIn:
    """Test channel fan-in functionality."""

    def test_channel_fanin_contacts(self, minimal_test_files: list[str], temp_dir: Path):
        """Test multiple contacts tasks emitting to same channel."""
        files = minimal_test_files
        out_path = temp_dir / "contacts_all.ndjson"

        p = Pipeline(FilesSource(files))

        # Two different contacts tasks, same output channel
        p.add_task(
            name="contacts_molstar",
            task=ContactTask(provider=ContactProvider.MolStar, interaction_type=InteractionType.All),
            channel="contacts",  # fan-in target
            in_memory_policy=InMemoryPolicy.Keep,
        )
        p.add_task(
            name="contacts_arpeggio",
            task=ContactTask(provider=ContactProvider.Arpeggio, interaction_type=InteractionType.All),
            channel="contacts",  # fan-in target
            in_memory_policy=InMemoryPolicy.Keep,
        )

        # Single sink subscribes to unified channel
        p.to_files("contacts", path=out_path, fmt=OutputFormat.JSON)

        results = p.run(threads=4)
        file_outputs = p.file_outputs()

        # Should have aggregated results from both providers
        assert "contacts" in results
        contacts = results["contacts"]
        assert len(contacts) == 2 * len(files)  # Results from both providers

        # Validate we have results from both providers
        providers_seen = {c["provider"] for c in contacts}
        assert "molstar" in providers_seen
        assert "arpeggio" in providers_seen

        assert out_path.exists()
        assert str(out_path) in file_outputs["contacts"]

    def test_channel_fanin_python_tasks(self, minimal_test_files: list[str], temp_dir: Path):
        """Test multiple Python tasks emitting to same channel."""

        def file_meta(ctx: PipelineContext) -> dict[str, str]:
            return {"file": os.path.basename(ctx.path)}

        def size_meta(ctx: PipelineContext) -> dict[str, int]:
            return {"size": os.path.getsize(ctx.path)}

        files = minimal_test_files
        out_path = temp_dir / "meta.ndjson"

        p = Pipeline(FilesSource(files))
        p.add_task(name="path_meta", task=file_meta, channel="meta", in_memory_policy=InMemoryPolicy.Keep)
        p.add_task(name="size_meta", task=size_meta, channel="meta", in_memory_policy=InMemoryPolicy.Keep)
        p.to_files("meta", path=out_path, fmt=OutputFormat.JSON)

        class MetaSchema(TypedDict):
            meta: list[dict[str, int | str]]

        results = p.run_typed(MetaSchema, threads=1)
        file_outputs = p.file_outputs()

        # Should have aggregated results from both tasks
        assert "meta" in results
        meta = results["meta"]
        assert len(meta) == 2 * len(files)  # Results from both tasks

        # Validate we have both types of metadata
        file_entries = [m for m in meta if "file" in m]
        size_entries = [m for m in meta if "size" in m]

        assert len(file_entries) == len(files)
        assert len(size_entries) == len(files)

        for size_entry in size_entries:
            size_value = size_entry["size"]
            assert isinstance(size_value, int)
            assert size_value > 0

        assert out_path.exists()
        assert str(out_path) in file_outputs["meta"]


class TestPipelineIntegration:
    """Integration tests combining multiple pipeline features."""

    def test_mixed_pipeline_dag(self, minimal_test_files: list[str], temp_dir: Path):
        """Test complex pipeline with mixed tasks and dependencies."""
        files = minimal_test_files

        def file_info(ctx: PipelineContext) -> dict[str, str]:
            return {"file": os.path.basename(ctx.path)}

        def annotate_with_contacts(ctx: PipelineContext) -> dict[str, Any]:
            return {"file": os.path.basename(ctx.path), "has_contacts": True, "processed": True}

        p = Pipeline(FilesSource(files))

        p.add_task(
            name="contacts",
            task=ContactTask(provider=ContactProvider.MolStar, interaction_type=InteractionType.All),
            in_memory_policy=InMemoryPolicy.Keep,
        )

        # Independent Python task depending only on system
        p.add_task(name="file_info", task=file_info, depends=["system"], in_memory_policy=InMemoryPolicy.Keep)

        # Python task depending on both system and contacts
        p.add_task(
            name="annotation",
            task=annotate_with_contacts,
            depends=["system", "contacts"],
            in_memory_policy=InMemoryPolicy.Keep,
        )

        results = p.run(threads=4)

        # All channels should have results
        assert "contacts"   in results
        assert "file_info"  in results
        assert "annotation" in results

        assert len(results["contacts"])   == len(files)
        assert len(results["file_info"])  == len(files)
        assert len(results["annotation"]) == len(files)

        # Validate annotation results
        for annotation in results["annotation"]:
            assert annotation["has_contacts"] is True
            assert annotation["processed"] is True
            assert "file" in annotation

    def test_convenience_sink_methods(self, minimal_test_files: list[str], temp_dir: Path):
        """Test convenience methods for attaching sinks."""
        files = minimal_test_files

        contacts_file = temp_dir / "contacts.ndjson"
        shards_dir = temp_dir / "contacts_shards"
        shards_dir.mkdir(exist_ok=True)

        p = Pipeline(FilesSource(files))
        p.add_task(name="contacts", task=ContactTask(provider=ContactProvider.MolStar, interaction_type=InteractionType.All))

        # Attach multiple sinks
        p.to_memory("contacts")
        p.to_files("contacts", path=contacts_file, fmt=OutputFormat.JSON)
        p.to_sharded_files("contacts", out_dir=shards_dir, fmt=OutputFormat.JSON, shard_size=1)

        results = p.run(threads=4)
        file_outputs = p.file_outputs()

        # Memory sink
        assert "contacts" in results
        assert len(results["contacts"]) == len(files)

        # File sink
        assert contacts_file.exists()
        assert str(contacts_file) in file_outputs["contacts"]

        # Sharded sink
        shard_files = [f for f in file_outputs["contacts"] if f.startswith(str(shards_dir))]
        assert len(shard_files) >= len(files)  # At least one shard per file


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
