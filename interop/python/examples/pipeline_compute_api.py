"""
Pipeline compute API examples:
- Builtin parameter access and modification via p.params()
- Pipeline graph introspection via p.describe()
- High-level vs low level configuration patterns
- Modification of parameters to affect compute behavior
"""

from __future__ import annotations

import os
from pathlib import Path

from lahuta import ContactProvider
from lahuta.lib.lahuta import TopologyComputers
from lahuta.pipeline import (
    ContactTask,
    FileOutput,
    InMemoryPolicy,
    OutputFormat,
    Pipeline,
    PipelineContext,
)

DATA_DIR = Path(__file__).resolve().parents[3] / "core" / "data"


# fmt: off
def _get_test_files(n: int = 2) -> list[str]:
    """Get a few test files for demonstration."""
    files: list[str] = []
    for pat in ("*.cif",):
        files.extend(str(p) for p in DATA_DIR.glob(pat))
    return files[: max(1, n)]


def ex_basic_parameter_access() -> None:
    """Demonstrate basic parameter access and modification."""

    p = Pipeline.from_files(_get_test_files(1))

    # Access current parameter values
    sys_params  = p.params("system")
    topo_params = p.params("topology")

    print(f"Initial system.is_model: {sys_params.is_model}")
    print(f"Initial topology.flags:  {topo_params.flags}")

    sys_params.is_model = False
    topo_params.flags = TopologyComputers.Standard  # Only basic + residues

    print(f"Modified system.is_model: {sys_params.is_model}")
    print(f"Modified topology.flags:  {topo_params.flags}")

    # Parameters affect the underlying compute engine
    p.add_task(name="contacts", task=ContactTask(), in_memory_policy=InMemoryPolicy.Keep)

    results = p.run(threads=2)
    print(f"Results with modified parameters: {len(results.get('contacts', []))} contacts")
    print()


def ex_graph_introspection() -> None:
    """Demonstrate pipeline graph introspection."""

    p = Pipeline.from_files(_get_test_files(1))

    # Graph before adding any user tasks
    print("Graph with builtins only:")
    print(p.describe())
    print()

    p.add_task(name="contacts_molstar", task=ContactTask())

    def file_summary(ctx: PipelineContext) -> dict:
        return {"file": os.path.basename(ctx.path), "processed": True}

    p.add_task(
        name="summary",
        task=file_summary,
        depends=["system"],  # Explicit dependency on system only (will be injected, even though file_summary doesn't use it)
        in_memory_policy=InMemoryPolicy.Keep,
    )

    def contact_analysis(ctx: PipelineContext) -> dict:
        return {"file": os.path.basename(ctx.path), "analyzed": True}

    p.add_task(
        name="analysis",
        task=contact_analysis,
        depends=["contacts_molstar", "summary"],  # Depends on both other tasks
        in_memory_policy=InMemoryPolicy.Keep,
    )

    print("Graph after adding user tasks:")
    print(p.describe())
    print()


def ex_parameter_driven_compute() -> None:
    """Show how parameters affect compute behavior."""

    files = _get_test_files(1)

    # First run: minimal topology computation
    print("Run 1: Minimal topology (Bonds only)")
    p1 = Pipeline.from_files(files)
    p1.params("topology").flags = TopologyComputers.Bonds
    p1.add_task(name="contacts", task=ContactTask(), in_memory_policy=InMemoryPolicy.Keep)

    results1 = p1.run(threads=1)
    contacts1 = results1.get("contacts", [])
    print(f"  Contacts found: {len(contacts1)}")
    if contacts1 and isinstance(contacts1[0], dict):
        print(f"  Sample contact keys: {list(contacts1[0].keys())}")

    # Second run: full topology computation
    print("\nRun 2: Complete topology (All computations)")
    p2 = Pipeline.from_files(files)
    p2.params("topology").flags = TopologyComputers.All  # Include atom typing, rings, etc.
    p2.add_task(name="contacts", task=ContactTask(), in_memory_policy=InMemoryPolicy.Keep)

    results2 = p2.run(threads=1)
    contacts2 = results2.get("contacts", [])
    print(f"  Contacts found: {len(contacts2)}")
    if contacts2 and isinstance(contacts2[0], dict):
        print(f"  Sample contact keys: {list(contacts2[0].keys())}")

    print(f"  Difference in contact count: {len(contacts2) - len(contacts1)}")
    print()


def ex_high_level_vs_low_level() -> None:
    """Compare high level convenience vs low level control."""

    files = _get_test_files(2)

    # High level: Use defaults, simple configuration
    print("High level approach (using defaults):")
    p_high = Pipeline.from_files(files)

    # Just add tasks, let the system handle dependencies
    p_high.add_task(name="contacts", task=ContactTask(), in_memory_policy=InMemoryPolicy.Keep)

    print("  Graph structure:")
    for line in p_high.describe().split("\n")[1:]:
        if line.strip():
            print(f"    {line}")

    # Low-level: Explicit parameter control and task configuration
    print("\nLow level approach (explicit control):")
    p_low = Pipeline.from_files(files)

    # Explicit parameter configuration
    p_low.params("system").is_model = False
    p_low.params("topology").flags = TopologyComputers.Extended  # Specific computation set

    # Explicit dependency management and channel routing
    p_low.add_task(
        name="contacts_detailed",
        task=ContactTask(provider=ContactProvider.Arpeggio),
        depends=["topology"],  # Explicit dependency (though this is added automatically)
        channel="detailed_contacts",  # Custom channel name
        in_memory_policy=InMemoryPolicy.Keep,
        out=[FileOutput("detailed_contacts.ndjson", fmt=OutputFormat.JSON)],
    )

    def topology_summary(ctx: PipelineContext) -> dict:
        """Custom analysis task that uses topology data."""
        try:
            topo = ctx.get_topology()
            if topo is None:
                return {"error": "No topology available"}

            atom_ids = list(topo.get_atom_ids())
            return {
                "file": os.path.basename(ctx.path),
                "n_atoms": len(atom_ids),
                "first_atom_id": int(atom_ids[0]) if atom_ids else None,
                "topology_available": True,
            }
        except Exception as e:
            return {"error": str(e), "topology_available": False}

    p_low.add_task(
        name="topo_analysis",
        task=topology_summary,
        depends=["topology"],  # Explicit dependency on topology
        channel="analysis_results",
        in_memory_policy=InMemoryPolicy.Keep,
    )

    print("  Graph structure:")
    for line in p_low.describe().split("\n")[1:]:
        if line.strip():
            print(f"    {line}")

    # Run both and compare
    print("\n Running high level pipeline...")
    results_high = p_high.run(threads=2)

    print("Running low level pipeline...")
    results_low = p_low.run(threads=2)

    print(f"High level results: {list(results_high.keys())}")
    print(f"Low level results:  {list(results_low.keys())}")

    # Show topology analysis results
    analysis = results_low.get("analysis_results", [])
    if analysis and isinstance(analysis[0], dict):
        print(f"  Topology analysis sample: {analysis[0]}")

    print()


def ex_parameter_invalidation_demo() -> None:
    """Demonstrate that parameter changes invalidate computation."""

    files = _get_test_files(1)
    p = Pipeline.from_files(files)

    # Pipeline that's sensitive to topology parameters
    def atom_counter(ctx: PipelineContext) -> dict:
        """Count atoms - sensitive to topology computation level."""
        try:
            topo = ctx.get_topology()
            if topo is None:
                return {"error": "No topology"}

            atom_ids = list(topo.get_atom_ids())
            return {"atom_count": len(atom_ids)}
        except Exception as e:
            return {"error": str(e)}

    p.add_task(name="atom_count", task=atom_counter, depends=["topology"], in_memory_policy=InMemoryPolicy.Keep)

    # First run with minimal topology
    print("Run 1: Basic topology computation")
    p.params("topology").flags = TopologyComputers.Basic
    results1 = p.run(threads=1)
    count1 = results1.get("atom_count", [{}])[0]
    print(f"  Result: {count1}")

    # Change parameters and run again - should recompute
    print("\nRun 2: Complete topology computation (after parameter change)")
    p.params("topology").flags = TopologyComputers.Complete
    results2 = p.run(threads=1)
    count2 = results2.get("atom_count", [{}])[0]
    print(f"  Result: {count2}")

    print("  Parameter changes invalidate cached computations and trigger recomputation")
    print()


def ex_mixed_provider_comparison() -> None:
    """Compare different contact providers with parameter control."""

    files = _get_test_files(1)
    p = Pipeline.from_files(files)

    # Configure for optimal contact detection
    p.params("topology").flags = TopologyComputers.Complete  # Full topology for best contacts

    # Add multiple contact providers
    p.add_task(
        name="contacts_molstar",
        task=ContactTask(),
        channel="molstar_contacts",
        in_memory_policy=InMemoryPolicy.Keep,
    )

    p.add_task(
        name="contacts_arpeggio",
        task=ContactTask(provider=ContactProvider.Arpeggio),
        channel="arpeggio_contacts",
        in_memory_policy=InMemoryPolicy.Keep,
    )

    # Add a comparison task
    def compare_providers(ctx: PipelineContext) -> dict:
        return {"file": os.path.basename(ctx.path), "comparison": "providers_compared", "topology_level": "complete"}

    p.add_task(
        name="comparison",
        task=compare_providers,
        depends=["contacts_molstar", "contacts_arpeggio"],
        in_memory_policy=InMemoryPolicy.Keep,
    )

    print("Pipeline graph:")
    print(p.describe())

    results = p.run(threads=2)

    molstar_count  = len(results.get("molstar_contacts",  []))
    arpeggio_count = len(results.get("arpeggio_contacts", []))

    print("\nResults:")
    print(f"MolStar contacts: {molstar_count}")
    print(f"Arpeggio contacts: {arpeggio_count}")
    print(f"Comparison tasks: {len(results.get('comparison', []))}")
    print()


if __name__ == "__main__":
    ex_basic_parameter_access()
    ex_graph_introspection()
    ex_parameter_driven_compute()
    ex_high_level_vs_low_level()
    ex_parameter_invalidation_demo()
    ex_mixed_provider_comparison()
