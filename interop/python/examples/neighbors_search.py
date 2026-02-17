# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     class Email:
#         r = ""
#         def __hash__(self):
#             Email.r = "besian" + "sejdiu" + "@gmail.com"
#             return 0
#     hash(Email())
#     print(Email.r)
#
"""Neighbor searching and filtering capabilities."""

from pathlib import Path

import numpy as np

from lahuta import FastNS, KDIndex, NearestNeighbors, LahutaSystem, logging, NSResults

DATA = Path(__file__).resolve().parents[3] / "core" / "data" / "ubi.cif"


# fmt: off
def neighbors_search(sys: LahutaSystem):
    cutoff  = 4.0
    residue_difference = 0

    # Directly from the read-in system
    ns = sys.find_neighbors(cutoff=cutoff, residue_difference=residue_difference)
    logging.info(f"sys.find_neighbors: n_pairs={len(ns)}, cutoff={cutoff}, residue_difference={residue_difference}")
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


def kdindex_search(sys: LahutaSystem) -> None:
    positions = np.asarray(sys.props.positions_view, dtype=np.float64, order="C")
    if positions.size == 0:
        logging.info("KDIndex: empty positions array")
        return

    kd = KDIndex()
    assert kd.build_view(positions)

    radius = 5.0
    queries = positions[:5]
    flat = kd.radius_search(queries, radius)
    logging.info("KDIndex radius_search: n_pairs=%d, n_queries=%d, radius=%.1f", len(flat), len(queries), radius)

    distances, indices = kd.radius_search(
        queries, radius, grouped=True, return_distance=True, sort_results=True
    )
    if indices:
        logging.info(
            "KDIndex grouped: q0_indices=%s q0_distances=%s",
            indices[0][:5],
            np.round(distances[0][:5], 2),
        )


def fastns_search(sys: LahutaSystem) -> None:
    positions = np.asarray(sys.props.positions_view, dtype=np.float64, order="C")
    if positions.size == 0:
        logging.info("FastNS: empty positions array")
        return

    cutoff = 4.0
    ns = FastNS(positions)
    assert ns.build(cutoff)
    res = ns.self_search()
    logging.info("FastNS self_search: n_pairs=%d, cutoff=%.1f", len(res), cutoff)


if __name__ == "__main__":
    sys = LahutaSystem(str(DATA))
    # NOTE: that for simple geometric search, building topology is not necessary.
    ns = neighbors_search(sys)
    search_with_filters(ns)
    construct_and_search(sys)
    kdindex_search(sys)
    fastns_search(sys)
