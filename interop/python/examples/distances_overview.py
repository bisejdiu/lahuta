"""Lahuta's distances and neighbor search APIs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from lahuta import LahutaSystem, NearestNeighbors, logging, metrics


@dataclass
class CoordinateSelection:
    """Holds a named coordinate selection takedn from a LahutaSystem."""

    name: str
    positions: np.ndarray  # shape (n, 3)


def load_ca_atoms(structure_path: Path) -> CoordinateSelection:
    system = LahutaSystem(str(structure_path))
    if not system.build_topology():
        raise RuntimeError(f"Failed to build topology for {structure_path}")

    props = system.props
    names = np.asarray([str(name).strip() for name in props.names])
    mask = names == "CA"
    ca_positions = props.positions[mask]
    if ca_positions.size == 0:
        raise RuntimeError("No Cα atoms found in the structure.")

    logging.info(
        "Loaded %d Cα atoms from %s (centroid=%s A)",
        ca_positions.shape[0],
        structure_path.name,
        np.round(ca_positions.mean(axis=0), 2),
    )
    return CoordinateSelection("Cα", ca_positions.astype(np.float64, copy=False))


def summarize_distribution(distances: np.ndarray, label: str) -> None:
    if distances.size == 0:
        logging.info("%s: no distances to summarise", label)
        return

    logging.info(
        "%s: min=%.2f A  median=%.2f A  mean=%.2f A  max=%.2f A\n",
        label,
        float(distances.min()),
        float(np.median(distances)),
        float(distances.mean()),
        float(distances.max()),
    )


def demo_pairwise_metrics(selection: CoordinateSelection) -> None:
    coords = selection.positions[:12]  # keep small for readability
    logging.info("Computing pairwise distances for %d %s atoms", coords.shape[0], selection.name)

    condensed = metrics.pdist(coords, squared=False)
    summarize_distribution(condensed, "metrics.pdist (Cα excerpt)")

    matrix = metrics.pairwise_distances(coords)
    logging.info("pairwise_distances matrix shape: %s", matrix.shape)
    logging.info("Row 0 distances (A): %s", np.round(matrix[0, 1:4], 2))


def demo_cross_distances(selection: CoordinateSelection) -> None:
    # Split the Cα cloud into two chunks and compute cross distances with cdist.
    coords = selection.positions
    mid = coords.shape[0] // 2
    ref = coords[:mid]
    qry = coords[mid : mid + 20]  # small probe set

    logging.info("Computing cross-distances between reference (%d) and query (%d) CA atoms", ref.shape[0], qry.shape[0])
    distances = metrics.cdist(qry, ref)
    logging.info("cdist result shape: %s", distances.shape)
    summarize_distribution(distances.ravel(), "metrics.cdist (CA vs CA)")


def demo_radius_neighbors(selection: CoordinateSelection, radius: float = 6.0) -> None:
    coords = selection.positions
    logging.info("Searching neighbours within %.1f A for %d CA atoms", radius, coords.shape[0])

    nn = NearestNeighbors(radius=radius, algorithm="kd_tree", sort_results=True).fit(coords)
    distances, indices = nn.radius_neighbors(return_distance=True)

    counts = np.array([len(idx) for idx in indices], dtype=np.int32)
    logging.info(
        "radius_neighbors: min=%d  median=%d  mean=%.1f  max=%d contacts / atom",
        int(counts.min()),
        int(np.median(counts)),
        float(counts.mean()),
        int(counts.max()),
    )

    if distances:
        summarize_distribution(np.concatenate(distances), "radius_neighbors distance distribution")

    logging.info(
        "Example neighbours for atom #0: indices=%s distances=%s A", indices[0][:5], np.round(distances[0][:5], 2)
    )


def select_structure_path() -> Path:
    """Locate the example structure shipped with the repository."""

    root = Path(__file__).resolve().parents[3]
    path = root / "data" / "ubi.cif"
    if not path.exists():
        raise FileNotFoundError(f"Example structure not found: {path}")
    return path


def main() -> None:
    logging.set_global_verbosity(logging.LogLevel.INFO)
    structure = select_structure_path()
    ca_selection = load_ca_atoms(structure)

    demo_pairwise_metrics(ca_selection)
    demo_cross_distances(ca_selection)
    demo_radius_neighbors(ca_selection)


if __name__ == "__main__":
    main()
