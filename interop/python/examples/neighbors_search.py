"""Neighbor searching and filtering capabilities."""

from pathlib import Path

import numpy as np

from lahuta import FastNS, LahutaSystem, NSResults, logging

DATA = Path(__file__).resolve().parents[3] / "data" / "ubi.cif"


# fmt: off
def neighbors_search(sys: LahutaSystem) -> NSResults:
    cutoff  = 4.0
    res_dif = 0

    ns = sys.find_neighbors(cutoff=cutoff, res_dif=res_dif).filter(cutoff)
    logging.info(f"neighbors: shape0={ns.pairs.shape[0]}, cutoff={cutoff}, res_dif={res_dif}")
    if not ns.pairs.shape[0]:
        pairs = ns.pairs[:5].tolist()
        logging.debug(f"pairs_head={pairs}")
    return ns


def search_with_filters(ns: NSResults) -> tuple[int, int, int]:
    ns_cut = ns.filter(2.5)  # Distance filter

    keep = list(range(10))
    ns_any = ns_cut.filter(keep)  # Keep by indices for either column

    ns_col0 = ns_cut.filter(keep, 0)  # filter first column, for indices in keep
    logging.info(f"neighbors: total={len(ns_cut)}, any_keep={len(ns_any)}, col0_keep={len(ns_col0)}")
    return len(ns_cut), len(ns_any), len(ns_col0)


def construct_and_search(sys: LahutaSystem) -> bool:
    ns = FastNS(sys.props.positions)
    ns.build(5)  # Build grid with a cutoff of 5.0 Angstroms
    res = ns.self_search()  # Self-search for neighbors
    logging.info(f"self_search: shape0={res.pairs.shape[0]}, cutoff=5.0")

    queries = np.array(
        [
            [0.2, 0.1, 0.0],
            [1.05, 0.05, 0.0],
        ],
        dtype=np.float64,
    )
    res = ns.search(queries)  # Search for neighbors of the queries
    logging.info(f"search: shape0={res.pairs.shape[0]}, queries={queries.shape}")

    return True


if __name__ == "__main__":
    sys = LahutaSystem(str(DATA))
    # NOTE: that for simple geometric search, building topology is not necessary.
    ns = neighbors_search(sys)
    search_with_filters(ns)
    construct_and_search(sys)
