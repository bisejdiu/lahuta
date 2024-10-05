from typing import TypeVar, Callable, Type, Any, Protocol, ClassVar

from lahuta._types.mdanalysis import AtomGroupType, UniverseType
from lahuta.core.topology.arc import Atom
from lahuta.core.topology.loaders import LahutaCPPLoader, TopologyLoader, Loader
from lahuta.lib import cLuni
from lahuta.lib._lahuta import IR, LahutaCPP
from lahuta.utils.radii import v_radii_assignment

# TODO:
# 1. Should _loaders store their __name__ or the actual class or implement a get_type method?
# 2. Write a factory to get the correct loader?


class _Loader(Protocol):
    def __init__(self, file_path: str) -> None:
        self.input: Loader

    def to_ir(self) -> IR: ...

    @staticmethod
    def from_ir(ir: IR) -> Loader: ...


class GemmiLoader:
    def __init__(self, file_path: str) -> None:  # -> LahutaCPP:
        print(f"Loading file: {file_path} using GemmiLoader")
        self.input = cLuni(file_path)

    def to_ir(self) -> IR:
        # Convert GemmiLoader to IR
        return IR(
            self.input.get_indices(),
            self.input.get_atomic_numbers(),
            self.input.get_names(),
            self.input.get_resids(),
            self.input.get_resnames(),
            self.input.get_chainlabels(),
            self.input.get_positions(),
        )

    @staticmethod
    def from_ir(ir: IR) -> LahutaCPP:
        # Convert IR to GemmiLoader
        return cLuni(ir)


class MDAnalysisLoader:
    def __init__(self, file_path: str) -> None:  # -> AtomGroupType:
        print(f"Loading file: {file_path} using MDAnalysisLoader")
        self.input = TopologyLoader(file_path).ag

    def to_ir(self) -> IR:
        # Convert MDAnalysisLoader to IR
        return IR(
            atom_indices=self.input.indices,
            atomic_numbers=self.input.atoms.atomicnums,
            atom_names=self.input.names,
            resids=self.input.resids,
            resnames=self.input.resnames,
            chainlabels=self.input.chainIDs.astype(str),
            positions=self.input.positions,
        )

    @staticmethod
    def from_ir(ir: IR) -> AtomGroupType:
        # Convert IR to MDAnalysisLoader
        import numpy as np
        import pandas as pd
        import MDAnalysis as mda

        chain_ids = pd.factorize(ir.chainlabels)[0]

        struct_arr = np.rec.fromarrays(  # type: ignore
            [ir.resnames, ir.resids, chain_ids],
            names=str("resnames, resids, chain_ids"),
        )

        resindices, uniques = pd.factorize(struct_arr)
        resnames, resids, chain_ids = uniques["resnames"], uniques["resids"], uniques["chain_ids"]

        uv: UniverseType = mda.Universe.empty(
            n_atoms=len(ir.atom_indices),
            n_residues=uniques.size,
            n_segments=len(ir.chainlabels),
            atom_resindex=resindices,
            residue_segindex=chain_ids,
            trajectory=True,
        )

        # Add topology attributes
        uv.add_TopologyAttr("names", ir.atom_names)
        # uv.add_TopologyAttr("type", self.arc.atoms.types)
        # uv.add_TopologyAttr("elements", ir.atom)
        # uv.add_TopologyAttr("vdw_radii", v_radii_assignment(self.arc.atoms.elements))
        uv.add_TopologyAttr("resnames", resnames)
        uv.add_TopologyAttr("resids", resids)
        # uv.add_TopologyAttr("chainIDs", self.arc.chains.auths)
        # uv.add_TopologyAttr("ids", self.arc.atoms.ids)
        # uv.add_TopologyAttr("tempfactors", self.arc.atoms.b_isos)

        # uv.atoms.positions = self.arc.atoms.coordinates
        # uv.filename = self.file_path

        return uv.atoms


# Loaders currently include: AtomGroupType, LahutaCPP


# Factory to get the correct loader
class LoaderFactory:
    loaders: ClassVar[dict[str, Type[_Loader]]] = {
        ".cif": GemmiLoader,
        ".pdb": MDAnalysisLoader,
    }
    # loaders: ClassVar[dict[str, Callable[[str], Loader]]] = {
    #     ".cif": GemmiLoader.load,
    #     ".pdb": MDAnalysisLoader.load,
    # }

    @staticmethod
    def load(file_path: str, extension: str) -> _Loader:
        if extension in LoaderFactory.loaders:
            return LoaderFactory.loaders[extension](file_path)

        raise ValueError(f"Unsupported file extension: {extension}")


# Conversion Logic
# S = TypeVar("S", bound=Loader)
# T = TypeVar("T", bound=IRData)
T = TypeVar("T", bound=_Loader)


class Converter:
    # converters: ClassVar[dict[tuple[Type[Any], Type[Any]], Callable[[Any], Any]]] = {}

    @staticmethod
    def convert(obj: _Loader, target_type: Type[T]) -> Loader:
        # Convert to intermediate representation
        if hasattr(obj, "to_ir") and hasattr(target_type, "from_ir"):
            ir = obj.to_ir()
            print(f"Converted to IR: {ir}")
            return target_type.from_ir(ir)

        raise ValueError(f"Conversion not supported: {type(obj).__name__} to {target_type.__name__}")


if __name__ == "__main__":
    file_name = "/Users/bsejdiu/projects/lahuta/cpp/data/1kx2_small.cif"

    # fmt: off
    # 1. Load a file of type .cif
    loader = LoaderFactory.load(file_name, ".cif")
    print(f"1.: {loader.input.names.shape}", loader.input)
    # 2. Convert the loaded file to IR
    ir_instance = loader.to_ir()
    print(f"2.: {len(ir_instance.atom_indices)}", ir_instance)
    # 3. Convert the IR instance to LahutaCPP
    lahuta_from_ir = loader.from_ir(ir_instance)
    print(f"3.: {lahuta_from_ir.names.shape}", lahuta_from_ir)
    # 4. Convert the LahutaCPP instance to AtomGroupType
    f2_instance = Converter.convert(loader, MDAnalysisLoader)
    _f2_instance = MDAnalysisLoader.from_ir(ir_instance)
    print(f"4.: {len(f2_instance.names)}", f2_instance)
    print(f"4.: {len(_f2_instance.names)}", _f2_instance)
    # 5. Convert AtomGroupType instance to IR
    # ir_instance = MDAnalysisLoader(file_name).to_ir()
    # print(f"5.: {len(ir_instance.atom_indices)}", ir_instance)
    # ir_instance = loader.to_ir()
    # 6. Convert IR instance to LahutaCPP
    lahuta_from_ir = loader.from_ir(ir_instance)
    print(f"6.: {lahuta_from_ir.names.shape}", lahuta_from_ir)

    # FIX: step 5 does not work because the current implementation does not allow for calling to_ir on the object created using to_ir
