"""
Demonstrates how to create a LahutaSystem from a file, build its topology,
and access various properties and perform neighbor searches.
"""

from pathlib import Path

from lahuta import ContactComputerType, LahutaSystem, TopologyBuildingOptions, logging


# fmt: off
def read_system(path: str | Path) -> LahutaSystem:
    sys = LahutaSystem(str(path))
    ok = sys.build_topology()
    if not ok:
        raise RuntimeError("Failed to build topology from the system.")
    return sys


def read_with_options(path: str | Path) -> LahutaSystem:
    logging.set_global_verbosity(logging.LogLevel.INFO)
    sys = LahutaSystem(str(path))
    opts = TopologyBuildingOptions()
    opts.cutoff = 4.5  # this is the cutoff for bond perception (4.5 is default)
    opts.compute_nonstandard_bonds = True
    opts.atom_typing_method = ContactComputerType.Molstar  # the default
    if not sys.build_topology(opts):
        raise RuntimeError("Failed to build topology with options from the system.")
    logging.info("Successfully built topology with options.")
    return sys


def show_properties(sys: LahutaSystem) -> None:
    p = sys.props
    logging.info(f"n_atoms  : {sys.n_atoms}")
    logging.info(f"file_name: {sys.file_name}")
    logging.info(f"positions.shape = {p.positions.shape} dtype={p.positions.dtype}")
    logging.info(f"indices.shape   = {p.indices.shape}   dtype={p.indices.dtype}")
    logging.info(f"names.shape     = {p.names.shape}     dtype={p.names.dtype}")

    logging.info(f"centroid [A]: {p.positions.mean(axis=0)}")

    backbone = {"N", "CA", "C", "O"}
    n_backbone = sum(1 for nm in p.names if str(nm).strip() in backbone)
    logging.info(f"backbone atoms: {n_backbone}\n")


def neighbor_search(sys: LahutaSystem, cutoff: float = 4.5, res_dif: int = 1) -> None:
    ns = sys.find_neighbors(cutoff=cutoff, res_dif=res_dif).filter(cutoff)

    dij = ns.get_sqrt_distances()

    if not dij.shape[0]:
        logging.info(f"No pairs within {cutoff:.2f} Å\n")
        return

    logging.info(f"pairs within {cutoff:.2f} Å: {dij.shape[0]}  (mean={float(dij.mean()):.3f}  min={float(dij.min()):.3f}  max={float(dij.max()):.3f})\n")


if __name__ == "__main__":
    data = Path(__file__).resolve().parents[3] / "data" / "ubi.cif"
    read_system(data)
    sys = read_with_options(data)
    show_properties(sys)
    neighbor_search(sys)
