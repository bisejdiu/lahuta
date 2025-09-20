"""Neighbor searching and filtering capabilities."""

from pathlib import Path

import numpy as np

from lahuta import NearestNeighbors, LahutaSystem, logging

DATA = Path(__file__).resolve().parents[3] / "data" / "ubi.cif"


# fmt: off
def neighbors_search(sys: LahutaSystem):
    cutoff  = 4.0
    res_dif = 0

    # Directly from the read-in system
    ns = sys.find_neighbors(cutoff=cutoff, res_dif=res_dif).filter(cutoff)
    logging.info(f"sys.find_neighbors: n_pairs={len(ns)}, cutoff={cutoff}, res_dif={res_dif}")
    return ns


def search_with_filters(ns: NSResults) -> tuple[int, int, int]:
    ns_cut = ns.filter(2.5)  # Distance filter

    keep = list(range(10))
    ns_any = ns_cut.filter(keep)  # Keep by indices for either column

    ns_col0 = ns_cut.filter(keep, 0)  # filter first column, for indices in keep
    logging.info(f"neighbors: total={len(ns_cut)}, any_keep={len(ns_any)}, col0_keep={len(ns_col0)}")
    return len(ns_cut), len(ns_any), len(ns_col0)


def construct_and_search(sys: LahutaSystem) -> bool:
    positions = sys.props.positions
    radius = 5.0
    nn = NearestNeighbors(radius=radius, algorithm="kd_tree").fit(positions)

    # Self neighbors
    idxs_self = nn.radius_neighbors(return_distance=False)
    n_pairs_est = (sum(map(len, idxs_self))) // 2
    logging.info(f"self radius_neighbors: approx unique_pairs={n_pairs_est}, radius={radius}")

    queries = np.array(
        [
            [0.2, 0.1, 0.0],
            [1.05, 0.05, 0.0],
        ],
        dtype=np.float64,
    )
    idxs = nn.radius_neighbors(queries, return_distance=False)
    logging.info(f"cross radius_neighbors: n_queries={len(idxs)}, radius={radius}")

    return True


if __name__ == "__main__":
    sys = LahutaSystem(str(DATA))
    # NOTE: that for simple geometric search, building topology is not necessary.
    ns = neighbors_search(sys)
    search_with_filters(ns)
    construct_and_search(sys)
