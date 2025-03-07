from typing import TypeVar, Callable, Type, Any, Protocol, ClassVar

from lahuta._types.mdanalysis import AtomGroupType, UniverseType
from lahuta.core.topology.loaders import LahutaCPPLoader, TopologyLoader, Loader
from lahuta.lib._lahuta import LahutaCPP as cLuni
from lahuta.lib._lahuta import factorize_residues
from lahuta.lib._lahuta import IR, LahutaCPP
from lahuta.utils.radii import v_radii_assignment
from lahuta.lib.utils import factorize

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
            self.input.positions,
        )

    @staticmethod
    def from_ir(ir: IR) -> LahutaCPP:
        # Convert IR to GemmiLoader
        return cLuni(ir)


class MDAnalysisLoader:
    def __init__(self, file_path: str) -> None:  # -> AtomGroupType:
        self.input = TopologyLoader(file_path).ag

    def to_ir(self) -> IR:
        # Convert MDAnalysisLoader to IR
        return IR(
            atom_indices=self.input.indices,
            atom_nums=self.input.atoms.atom_nums,
            atom_names=self.input.names,
            resids=self.input.resids,
            resnames=self.input.resnames,
            chainlabels=self.input.chainIDs.astype(str),
            positions=self.input.positions,
        )

    @staticmethod
    def from_ir(ir: IR) -> AtomGroupType:
        # Convert IR to MDAnalysisLoader
        import MDAnalysis as mda

        unique_residues = factorize(ir.resnames, ir.resids, ir.chainlabels)
        uv: UniverseType = mda.Universe.empty(
            n_atoms=len(ir.atom_indices),
            n_residues=len(unique_residues.chains),
            n_segments=cLuni.count_unique(ir.chainlabels),
            atom_resindex=unique_residues.resindices,
            residue_segindex=factorize(unique_residues.chains),
            # residue_segindex=cLuni.factorize(rrc.chains),
            trajectory=True,
        )

        # Add topology attributes
        elements = cLuni.find_elements(ir.atom_nums)
        uv.add_TopologyAttr("names", ir.atom_names)
        # uv.add_TopologyAttr("type", self.arc.atoms.types)
        # uv.add_TopologyAttr("elements", ir.atom)
        uv.add_TopologyAttr("elements", elements)
        uv.add_TopologyAttr("resnames", unique_residues.resnames)
        uv.add_TopologyAttr("resids", unique_residues.resids)
        uv.add_TopologyAttr("chainIDs", ir.chainlabels)
        uv.add_TopologyAttr("ids", ir.atom_indices)

        # Make faster?
        uv.add_TopologyAttr("vdw_radii", v_radii_assignment(elements))

        # FIX: we need to add tempfactors
        # uv.add_TopologyAttr("tempfactors", self.arc.atoms.b_isos)

        uv.atoms.positions = ir.positions

        return uv.atoms


class LoaderFactory:
    loaders: ClassVar[dict[str, Type[_Loader]]] = {
        "cif": GemmiLoader,
        "cif.gz": GemmiLoader,
        "pdb": MDAnalysisLoader,
    }

    @staticmethod
    def load(file_path: str, extension: str) -> _Loader:
        if extension in LoaderFactory.loaders:
            return LoaderFactory.loaders[extension](file_path)

        raise ValueError(f"Unsupported file extension: {extension}")


T = TypeVar("T", bound=_Loader)


class Converter:
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
    file_name = "/Users/bsejdiu/projects/lahuta/cpp/build/mod_1kx2.cif"

    # fmt: off
    # 1. Load a file of type .cif
    import numpy as np
    loader = LoaderFactory.load(file_name, ".cif")
    print("xx: ", np.unique(loader.input.chainlabels))
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
