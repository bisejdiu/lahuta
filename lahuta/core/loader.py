from typing import TypeVar, Callable, Type, Any, Protocol, ClassVar

from lahuta._types.mdanalysis import AtomGroupType
from lahuta.core.topology.arc import Atom
from lahuta.core.topology.loaders import LahutaCPPLoader, TopologyLoader, Loader
from lahuta.lib import cLuni
from lahuta.lib._lahuta import IR, LahutaCPP


# class IRData(Protocol):
#     def __init__(self):
#         self.input: LahutaCPP
#
#     # @staticmethod
#     def to_ir() -> IR: ...
#
#     @staticmethod
#     def from_ir(ir: IR) -> _Loader: ...


class _Loader(Protocol):
    def __init__(self, file_path: str) -> None:
        self.input: Loader

    # @staticmethod
    # def load(file_path: str) -> Any: ...

    # @staticmethod
    def to_ir(self) -> IR: ...

    @staticmethod
    def from_ir(ir: IR) -> Loader: ...


class GemmiLoader:
    def __init__(self, file_path: str) -> None:  # -> LahutaCPP:
        print(f"Loading file: {file_path} using GemmiLoader")
        self.input = cLuni(file_path)
        # return cLuni(file_path)

    # @staticmethod
    def to_ir(self) -> IR:
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
        # Convert IRType data to GemmiLoader
        ...


class MDAnalysisLoader:
    # @staticmethod
    def __init__(self, file_path: str) -> None:  # -> AtomGroupType:
        print(f"Loading file: {file_path} using MDAnalysisLoader")
        self.input = TopologyLoader(file_path).ag
        # return TopologyLoader(file_path).ag

    # @staticmethod
    def to_ir(self) -> IR:
        # Convert MDAnalysisLoader data to IR
        ...
        # return IR(
        #     atom_indices=ag.indices,
        #     atomic_numbers=ag.atoms.atomic_numbers,
        #     atom_names=ag.names,
        #     resids=ag.resids,
        #     resnames=ag.resnames,
        #     chainlabels=ag.chainlabels,
        #     positions=ag.positions,
        # )

    @staticmethod
    def from_ir(ir: IR) -> AtomGroupType:
        # Convert IRType data to MDAnalysisLoader
        ...


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

    # @staticmethod
    # def register_conversion(source: Type[S], target: Type[T], func: Callable[[S], T]) -> None:
    #     Converter.converters[(source, target)] = func

    @staticmethod
    def convert(obj: _Loader, target_type: Type[T]) -> T:
        # Convert to intermediate representation
        if hasattr(obj, "to_ir") and hasattr(target_type, "from_ir"):
            ir = obj.to_ir()
            print(f"Converted to IR: {ir}: {ir.atom_indices}")
            return target_type.from_ir(ir)

        raise ValueError(f"Conversion not supported: {type(obj).__name__} to {target_type.__name__}")

    # @staticmethod
    # def convert(obj: Loader, target_type: Type[T]) -> T:
    #     source_type = type(obj)
    #     func = Converter.converters.get((source_type, target_type))
    #     if func:
    #         return func(obj)
    #     raise ValueError(f"No conversion registered from {source_type.__name__} to {target_type.__name__}")


# Register conversions
# def f2_to_f1(f2: AtomGroupType) -> LahutaCPP:
#     print("Converting AtomGroupType to LahutaCPP")
#     return cLuni(f2.file_name)
#
#
# def f1_to_f2(f1: LahutaCPP) -> AtomGroupType:
#     print("Converting LahutaCPP to AtomGroupType")
#     return
#     return TopologyLoader(f1.file_name).ag


# Converter.register_conversion(AtomGroupType, LahutaCPP, f2_to_f1)
# Converter.register_conversion(LahutaCPP, AtomGroupType, f1_to_f2)

# Example Usage
if __name__ == "__main__":
    # Load a file of type .ext1
    # f1_instance = LoaderFactory.load("f1.ext1", ".ext1")
    file_name = "/Users/bsejdiu/projects/lahuta/cpp/data/1kx2_small.cif"
    f1_instance = LoaderFactory.load(file_name, ".cif")
    print(f"F1 compute result: {f1_instance.input.names}")

    # Convert LahutaCPP instance to AtomGroupType
    f1_instance = Converter.convert(f1_instance, MDAnalysisLoader)
    # print(f"F1 compute result: {f1_instance.input.names}")

    # Convert F1Type instance to F2Type
    # f2_instance = Converter.convert(f1_instance, AtomGroupType)
    # print(f"F2 compute result: {f2_instance.names}")
