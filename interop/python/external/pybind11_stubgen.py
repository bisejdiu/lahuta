"""
Copyright (c) 2019 Sergei Izmailov <sergei.a.izmailov@gmail.com>, All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to the author of this software, without
imposing a separate written license agreement for such Enhancements, then you
hereby grant the following license: a non-exclusive, royalty-free perpetual
license to install, use, modify, prepare derivative works, incorporate into
other computer software, distribute, and sublicense such enhancements or
derivative works thereof, in binary and source code form.
"""

from __future__ import annotations

import abc
import ast
import builtins
import dataclasses
import importlib
import inspect
import logging
import re
import sys
import types
from argparse import ArgumentParser, Namespace
from collections.abc import Sequence
from dataclasses import dataclass
from dataclasses import field as field_
from logging import getLogger
from pathlib import Path
from typing import Any, Tuple, Union


if sys.version_info >= (3, 8):
    from typing import Literal

    Modifier = Literal["static", "class", None]
else:
    from typing import Optional

    Modifier = Optional[str]


class Identifier(str):
    pass


class Decorator(str):
    pass


class Docstring(str):
    pass


@dataclass
class InvalidExpression:
    text: str

    def __str__(self):
        return f"Invalid python expression `{self.text}`"


class QualifiedName(Tuple[Identifier, ...]):
    @classmethod
    def from_str(cls, name: str) -> QualifiedName:
        return QualifiedName(Identifier(part) for part in name.split("."))

    def __str__(self):
        return ".".join(self)

    @property
    def parent(self) -> QualifiedName:
        return QualifiedName(self[:-1])


@dataclass
class Value:
    repr: str
    is_print_safe: bool = False

    def __str__(self):
        return self.repr


@dataclass
class ResolvedType:
    name: QualifiedName
    parameters: list[ResolvedType | Value | InvalidExpression] | None = field_(default=None)

    def __str__(self):
        if self.parameters:
            param_str = "[" + ", ".join(str(p) for p in self.parameters) + "]"
        else:
            param_str = ""
        return f"{self.name}{param_str}"


@dataclass
class Alias:
    name: Identifier
    origin: QualifiedName


Annotation = Union[ResolvedType, Value, InvalidExpression]


@dataclass
class TypeVar_:
    name: Identifier
    constraints: list[Annotation] = field_(default_factory=list)
    bound: Annotation | None = field_(default=None)
    covariant: bool = field_(default=False)
    contravariant: bool = field_(default=False)

    def __str__(self) -> str:
        return (
            f'{self.name} = typing.TypeVar("{self.name}"'
            + (", " if self.constraints else "")
            + (", ".join(str(c) for c in self.constraints))
            + (f", bound={str(self.bound)}" if self.bound is not None else "")
            + (f", covariant={self.covariant}" if self.covariant else "")
            + (f", contravariant={self.contravariant}" if self.contravariant else "")
            + ")"
        )


@dataclass
class Attribute:
    name: Identifier
    value: Value | None
    annotation: Annotation | None = field_(default=None)


@dataclass
class Argument:
    name: Identifier | None
    pos_only: bool = field_(default=False)
    kw_only: bool = field_(default=False)
    variadic: bool = field_(default=False)
    kw_variadic: bool = field_(default=False)
    default: Value | InvalidExpression | None = field_(default=None)
    annotation: Annotation | None = field_(default=None)

    def __str__(self):
        result = []
        if self.variadic:
            result += ["*"]
        if self.kw_variadic:
            result += ["**"]
        result += [f"{self.name}"]
        if self.annotation:
            result += [f": {self.annotation}"]
        if self.default:
            result += [f" = {self.default}"]

        return "".join(result)


@dataclass
class Function:
    name: Identifier
    args: list[Argument] = field_(default_factory=list)
    returns: Annotation | None = field_(default=None)
    doc: Docstring | None = field_(default=None)
    decorators: list[Decorator] = field_(default_factory=list)

    def __str__(self):
        return f"{self.name}({', '.join(str(arg) for arg in self.args)}) -> {self.returns}"


@dataclass
class Property:
    name: Identifier
    modifier: Modifier
    doc: Docstring | None = field_(default=None)
    getter: Function | None = field_(default=None)
    setter: Function | None = field_(default=None)


@dataclass
class Method:
    function: Function
    modifier: Modifier


@dataclass
class Field:
    attribute: Attribute
    modifier: Modifier


@dataclass
class Class:
    name: Identifier
    doc: Docstring | None = field_(default=None)
    bases: list[QualifiedName] = field_(default_factory=list)
    classes: list[Class] = field_(default_factory=list)
    fields: list[Field] = field_(default_factory=list)
    methods: list[Method] = field_(default_factory=list)
    properties: list[Property] = field_(default_factory=list)
    aliases: list[Alias] = field_(default_factory=list)


@dataclass(eq=True, frozen=True)
class Import:
    name: Identifier | None
    origin: QualifiedName


@dataclass
class Module:
    name: Identifier
    doc: Docstring | None = field_(default=None)
    classes: list[Class] = field_(default_factory=list)
    functions: list[Function] = field_(default_factory=list)
    sub_modules: list[Module] = field_(default_factory=list)
    attributes: list[Attribute] = field_(default_factory=list)
    imports: set[Import] = field_(default_factory=set)
    aliases: list[Alias] = field_(default_factory=list)
    type_vars: list[TypeVar_] = field_(default_factory=list)
    is_package: bool = field_(default=False)


class FixedSize:
    def __init__(self, *dim: int):
        self.dim: tuple[int, ...] = dim

    def __repr__(self):
        return f"{self.__module__}.{self.__class__.__qualname__}({', '.join(str(d) for d in self.dim)})"


class DynamicSize:
    def __init__(self, *dim: int | str):
        self.dim: tuple[int | str, ...] = dim

    def __repr__(self):
        return f"{self.__module__}.{self.__class__.__qualname__}({', '.join(repr(d) for d in self.dim)})"


class ParserError(Exception):
    pass


class InvalidIdentifierError(ParserError):
    def __init__(self, name: Identifier, found_at: QualifiedName):
        super().__init__()
        self.name = name
        self.at = found_at

    def __str__(self):
        return f"Invalid identifier '{self.name}' at '{self.at}'"


class InvalidExpressionError(ParserError):
    def __init__(self, expression: str):
        super().__init__()
        self.expression: str = expression

    def __str__(self):
        return f"Invalid expression '{self.expression}'"


class NameResolutionError(ParserError):
    def __init__(self, name: QualifiedName):
        super().__init__()
        self.name = name

    def __str__(self):
        return f"Can't find/import '{self.name}'"


class IParser(abc.ABC):
    @abc.abstractmethod
    def handle_alias(self, path: QualifiedName, origin: Any) -> Alias | None: ...

    @abc.abstractmethod
    def handle_attribute(self, path: QualifiedName, attr: Any) -> Attribute | None: ...

    @abc.abstractmethod
    def handle_bases(self, path: QualifiedName, bases: tuple[type, ...]) -> list[QualifiedName]: ...

    @abc.abstractmethod
    def handle_class(self, path: QualifiedName, class_: type) -> Class | None: ...

    @abc.abstractmethod
    def handle_class_member(
        self, path: QualifiedName, class_: type, obj: Any
    ) -> Docstring | Alias | Class | list[Method] | Field | Property | None: ...

    @abc.abstractmethod
    def handle_docstring(self, path: QualifiedName, doc: Any) -> Docstring | None: ...

    @abc.abstractmethod
    def handle_field(self, path: QualifiedName, field: Any) -> Field | None: ...

    @abc.abstractmethod
    def handle_function(self, path: QualifiedName, func: Any) -> list[Function]: ...

    @abc.abstractmethod
    def handle_import(self, path: QualifiedName, origin: Any) -> Import | None: ...

    @abc.abstractmethod
    def handle_method(self, path: QualifiedName, method: Any) -> list[Method]: ...

    @abc.abstractmethod
    def handle_module(self, path: QualifiedName, module: types.ModuleType) -> Module | None: ...

    @abc.abstractmethod
    def handle_module_member(
        self, path: QualifiedName, module: types.ModuleType, obj: Any
    ) -> Docstring | Import | Alias | Class | list[Function] | Attribute | Module | None: ...

    @abc.abstractmethod
    def handle_property(self, path: QualifiedName, prop: Any) -> Property | None: ...

    @abc.abstractmethod
    def handle_type(self, type_: type) -> QualifiedName: ...

    @abc.abstractmethod
    def handle_value(self, value: Any) -> Value: ...

    @abc.abstractmethod
    def parse_args_str(self, args_str: str) -> list[Argument]: ...

    @abc.abstractmethod
    def parse_annotation_str(self, annotation_str: str) -> ResolvedType | InvalidExpression | Value: ...

    @abc.abstractmethod
    def parse_value_str(self, value: str) -> Value | InvalidExpression: ...

    @abc.abstractmethod
    def report_error(self, error: ParserError) -> None: ...

    @abc.abstractmethod
    def finalize(self): ...


class LocalErrors:
    def __init__(self, path: QualifiedName, errors: set[str], stack: list[LocalErrors]):
        self.path: QualifiedName = path
        self.errors: set[str] = errors
        self._stack: list[LocalErrors] = stack

    def __enter__(self) -> LocalErrors:
        self._stack.append(self)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        top = self._stack.pop()
        assert top == self


class LoggerData(IParser):
    def __init__(self):
        super().__init__()
        self.stack: list[LocalErrors] = []

    def __new_layer(self, path: QualifiedName) -> LocalErrors:
        return LocalErrors(path, errors=set(), stack=self.stack)

    def handle_module(self, path: QualifiedName, module: types.ModuleType) -> Module | None:
        with self.__new_layer(path):
            return super().handle_module(path, module)

    def handle_class(self, path: QualifiedName, class_: type) -> Class | None:
        with self.__new_layer(path):
            return super().handle_class(path, class_)

    def handle_function(self, path: QualifiedName, class_: type) -> list[Function]:
        with self.__new_layer(path):
            return super().handle_function(path, class_)

    def handle_method(self, path: QualifiedName, class_: type) -> list[Method]:
        with self.__new_layer(path):
            return super().handle_method(path, class_)

    @property
    def current_path(self) -> QualifiedName:
        assert len(self.stack) != 0
        return self.stack[-1].path

    @property
    def reported_errors(self) -> set[str]:
        assert len(self.stack) != 0
        return self.stack[-1].errors


logger = getLogger("pybind11_stubgen")


class LogErrors(IParser):
    def report_error(self: LoggerData, error: ParserError) -> None:
        error_str = f"In {self.current_path} : {error}"
        if error_str not in self.reported_errors:
            logger.error(error_str)
            self.reported_errors.add(error_str)
        super().report_error(error)


class IgnoreUnresolvedNameErrors(IParser):
    def __init__(self):
        super().__init__()
        self.__regex: re.Pattern | None = None

    def set_ignored_unresolved_names(self, regex: re.Pattern):
        self.__regex = regex

    def report_error(self, error: ParserError):
        if self.__regex is not None:
            if isinstance(error, NameResolutionError):
                if self.__regex.match(str(error.name)):
                    return
        super().report_error(error)


class IgnoreInvalidExpressionErrors(IParser):
    def __init__(self):
        super().__init__()
        self.__regex: re.Pattern | None = None

    def set_ignored_invalid_expressions(self, regex: re.Pattern):
        self.__regex = regex

    def report_error(self, error: ParserError):
        if self.__regex is not None:
            if isinstance(error, InvalidExpressionError):
                if self.__regex.match(error.expression):
                    return
        super().report_error(error)


class IgnoreInvalidIdentifierErrors(IParser):
    def __init__(self):
        super().__init__()
        self.__regex: re.Pattern | None = None

    def set_ignored_invalid_identifiers(self, regex: re.Pattern):
        self.__regex = regex

    def report_error(self, error: ParserError):
        if self.__regex is not None:
            if isinstance(error, InvalidIdentifierError):
                if self.__regex.match(str(error.name)):
                    return
        super().report_error(error)


class IgnoreAllErrors(IParser):
    def report_error(self, error: ParserError):
        return None


class SuggestCxxSignatureFix(IParser):
    def __init__(self):
        super().__init__()
        self.suggest_cxx_fix = False

    def report_error(self, error: ParserError):
        if isinstance(error, InvalidExpressionError):
            expression = error.expression
            if "::" in expression or expression.endswith(">"):
                self.suggest_cxx_fix = True

        super().report_error(error)

    def finalize(self):
        if self.suggest_cxx_fix:
            logger.warning(
                "Raw C++ types/values were found in signatures extracted "
                "from docstrings.\n"
                "Please check the corresponding sections of pybind11 documentation "
                "to avoid common mistakes in binding code:\n"
                " - https://pybind11.readthedocs.io/en/latest/advanced/misc.html"
                "#avoiding-cpp-types-in-docstrings\n"
                " - https://pybind11.readthedocs.io/en/latest/advanced/functions.html"
                "#default-arguments-revisited"
            )
        super().finalize()


class TerminateOnFatalErrors(IParser):
    def __init__(self):
        super().__init__()
        self.__found_fatal_errors = False

    def report_error(self, error: ParserError):
        self.__found_fatal_errors = True
        super().report_error(error)

    def finalize(self):
        super().finalize()
        if self.__found_fatal_errors:
            logger.info("Terminating due to previous errors")
            sys.exit(1)


class FilterTypingModuleAttributes(IParser):
    __sentinel = object()

    def handle_attribute(self, path: QualifiedName, attr: Any) -> Attribute | None:
        import typing

        if getattr(typing, str(path[-1]), self.__sentinel) is attr:
            return None
        return super().handle_attribute(path, attr)


class FilterClassMembers(IParser):
    __attribute_blacklist: set[Identifier] = {
        *map(
            Identifier,
            (
                "__annotations__",
                "__builtins__",
                "__cached__",
                "__file__",
                "__firstlineno__",
                "__loader__",
                "__name__",
                "__package__",
                "__path__",
                "__spec__",
                "__static_attributes__",
            ),
        )
    }
    __class_member_blacklist: set[Identifier] = {
        *map(
            Identifier,
            (
                "__annotations__",
                "__class__",
                "__dict__",
                "__module__",
                "__qualname__",
                "__weakref__",
            ),
        )
    }

    __method_blacklist: set[Identifier] = {*map(Identifier, ("__dir__", "__sizeof__"))}

    def handle_class_member(
        self, path: QualifiedName, class_: type, member: Any
    ) -> Docstring | Alias | Class | list[Method] | Field | Property | None:
        name = path[-1]
        if name in self.__class_member_blacklist:
            return None
        if not hasattr(class_, "__dict__") or name not in class_.__dict__:
            return None
        return super().handle_class_member(path, class_, member)

    def handle_method(self, path: QualifiedName, value: Any) -> list[Method]:
        if path[-1] in self.__method_blacklist:
            return []
        return super().handle_method(path, value)

    def handle_attribute(self, path: QualifiedName, attr: Any) -> Attribute | None:
        if path[-1] in self.__attribute_blacklist:
            return None
        return super().handle_attribute(path, attr)


class FilterPybindInternals(IParser):
    __attribute_blacklist: set[Identifier] = {*map(Identifier, ("__entries",))}

    __class_blacklist: set[Identifier] = {*map(Identifier, ("pybind11_type",))}

    def handle_attribute(self, path: QualifiedName, value: Any) -> Attribute | None:
        if path[-1] in self.__attribute_blacklist:
            return None
        return super().handle_attribute(path, value)

    def handle_class_member(
        self, path: QualifiedName, class_: type, member: Any
    ) -> Docstring | Alias | Class | list[Method] | Field | Property | None:
        name = path[-1]
        if name in self.__class_blacklist:
            return None
        if name.startswith("__pybind11_module"):
            return None
        if name.startswith("_pybind11_conduit_v1_"):
            return None
        return super().handle_class_member(path, class_, member)


class FilterInvalidIdentifiers(IParser):
    def handle_module_member(
        self, path: QualifiedName, module: types.ModuleType, obj: Any
    ) -> Docstring | Import | Alias | Class | list[Function] | Attribute | Module | None:
        if not path[-1].isidentifier():
            self.report_error(InvalidIdentifierError(path[-1], path.parent))
            return None
        return super().handle_module_member(path, module, obj)

    def handle_class_member(
        self, path: QualifiedName, class_: type, obj: Any
    ) -> Docstring | Alias | Class | list[Method] | Field | Property | None:
        if not path[-1].isidentifier():
            self.report_error(InvalidIdentifierError(path[-1], path.parent))
            return None
        return super().handle_class_member(path, class_, obj)


class FilterPybind11ViewClasses(IParser):
    def handle_module_member(
        self, path: QualifiedName, module: types.ModuleType, obj: Any
    ) -> Docstring | Import | Alias | Class | list[Function] | Attribute | Module | None:
        result = super().handle_module_member(path, module, obj)

        if isinstance(result, Class) and str(result.name) in [
            "ItemsView",
            "KeysView",
            "ValuesView",
        ]:
            return None

        return result


logger2 = getLogger("pybind11_stubgen")


class RemoveSelfAnnotation(IParser):
    __any_t_name = QualifiedName.from_str("Any")
    __typing_any_t_name = QualifiedName.from_str("typing.Any")

    def handle_method(self, path: QualifiedName, method: Any) -> list[Method]:
        methods = super().handle_method(path, method)
        for method in methods:
            self._remove_self_arg_annotation(path, method.function)
        return methods

    def handle_property(self, path: QualifiedName, prop: Any) -> Property | None:
        prop = super().handle_property(path, prop)
        if prop is not None:
            if prop.getter is not None:
                self._remove_self_arg_annotation(path, prop.getter)
            if prop.setter is not None:
                self._remove_self_arg_annotation(path, prop.setter)

        return prop

    def _remove_self_arg_annotation(self, path: QualifiedName, func: Function) -> None:
        if len(func.args) == 0:
            return
        fully_qualified_class_name = QualifiedName(path[:-1])
        first_arg = func.args[0]
        if (
            first_arg.name == "self"
            and isinstance(first_arg.annotation, ResolvedType)
            and not first_arg.annotation.parameters
            and (
                first_arg.annotation.name
                in {
                    self.__any_t_name,
                    self.__typing_any_t_name,
                    fully_qualified_class_name,
                    fully_qualified_class_name[-len(first_arg.annotation.name) :],
                }
            )
        ):
            first_arg.annotation = None


class FixMissingImports(IParser):
    def __init__(self):
        super().__init__()
        self.__extra_imports: set[Import] = set()
        self.__current_module: types.ModuleType | None = None
        self.__current_class: type | None = None

    def handle_alias(self, path: QualifiedName, origin: Any) -> Alias | None:
        result = super().handle_alias(path, origin)
        if result is None:
            return None
        self._add_import(result.origin)
        return result

    def handle_attribute(self, path: QualifiedName, attr: Any) -> Attribute | None:
        result = super().handle_attribute(path, attr)
        if result is None:
            return None
        if isinstance(result.annotation, ResolvedType):
            self._add_import(result.annotation.name)
        return result

    def handle_class(self, path: QualifiedName, class_: type) -> Class | None:
        old_class = self.__current_class
        self.__current_class = class_
        result = super().handle_class(path, class_)
        self.__current_class = old_class
        return result

    def handle_import(self, path: QualifiedName, origin: Any) -> Import | None:
        result = super().handle_import(path, origin)
        if result is None:
            return None
        self.__extra_imports.add(result)
        return result

    def handle_module(self, path: QualifiedName, module: types.ModuleType) -> Module | None:
        old_imports = self.__extra_imports
        old_module = self.__current_module
        self.__extra_imports = set()
        self.__current_module = module
        result = super().handle_module(path, module)
        if result is not None:
            result.imports |= self.__extra_imports
        self.__extra_imports = old_imports
        self.__current_module = old_module
        return result

    def handle_type(self, type_: type) -> QualifiedName:
        result = super().handle_type(type_)
        if not inspect.ismodule(type_):
            self._add_import(result)
        return result

    def handle_value(self, value: Any) -> Value:
        result = super().handle_value(value)
        if inspect.isroutine(value) and result.is_print_safe:
            self._add_import(QualifiedName.from_str(result.repr))
        return result

    def parse_annotation_str(self, annotation_str: str) -> ResolvedType | InvalidExpression | Value:
        result = super().parse_annotation_str(annotation_str)
        if isinstance(result, ResolvedType):
            self._add_import(result.name)
        return result

    def _add_import(self, name: QualifiedName) -> None:
        if len(name) == 0:
            return
        if len(name) == 1 and len(name[0]) == 0:
            return
        if hasattr(builtins, name[0]):
            return
        if self.__current_class is not None and hasattr(self.__current_class, name[0]):
            return
        if self.__current_module is not None and hasattr(self.__current_module, name[0]):
            return
        module_name = self._get_parent_module(name)
        if module_name is None:
            self.report_error(NameResolutionError(name))
            return
        self.__extra_imports.add(Import(name=None, origin=module_name))

    def _get_parent_module(self, name: QualifiedName) -> QualifiedName | None:
        parent = name.parent
        while len(parent) != 0:
            if self._is_module(parent):
                if not self._is_accessible(name, from_module=parent):
                    return None
                return parent
            parent = parent.parent
        return None

    def _is_module(self, name: QualifiedName):
        try:
            return importlib.import_module(str(name)) is not None
        except ModuleNotFoundError:
            return False

    def _is_accessible(self, name: QualifiedName, from_module: QualifiedName) -> bool:
        try:
            parent = importlib.import_module(str(from_module))
        except ModuleNotFoundError:
            return False
        relative_path = name[len(from_module) :]
        for part in relative_path:
            if not hasattr(parent, part):
                return False
            parent = getattr(parent, part)
        return True


class FixMissing__future__AnnotationsImport(IParser):
    def handle_module(self, path: QualifiedName, module: types.ModuleType) -> Module | None:
        result = super().handle_module(path, module)
        if result is None:
            return None
        result.imports.add(self._future(Identifier("annotations")))
        return result

    def _future(self, feature: Identifier) -> Import:
        return Import(name=feature, origin=QualifiedName((Identifier("__future__"), feature)))


class FixMissing__all__Attribute(IParser):
    def handle_module(self, path: QualifiedName, module: types.ModuleType) -> Module | None:
        result = super().handle_module(path, module)
        if result is None:
            return None

        for attr in result.attributes:
            if attr.name == Identifier("__all__"):
                return result

        all_names: list[str] = sorted(
            set(
                filter(
                    lambda name: not name.startswith("_"),
                    map(
                        str,
                        (
                            *(class_.name for class_ in result.classes),
                            *(attr.name for attr in result.attributes),
                            *(func.name for func in result.functions),
                            *(alias.name for alias in result.aliases),
                            *(import_.name for import_ in result.imports if import_.name is not None),
                            *(sub_module.name for sub_module in result.sub_modules),
                        ),
                    ),
                )
            )
        )

        result.attributes.append(
            Attribute(
                name=Identifier("__all__"),
                value=self.handle_value(all_names),
                annotation=ResolvedType(name=QualifiedName.from_str("list[str]")),
            )
        )

        return result


class FixBuiltinTypes(IParser):
    _any_type = QualifiedName.from_str("typing.Any")

    def handle_type(self, type_: type) -> QualifiedName:
        if type_.__qualname__ == "PyCapsule" and type_.__module__ == "builtins":
            return self._any_type

        result = super().handle_type(type_)

        if result[0] == "builtins":
            if result[1] == "NoneType":
                return QualifiedName((Identifier("None"),))
            if result[1] in ("function", "builtin_function_or_method"):
                callable_t = self.parse_annotation_str("typing.Callable")
                assert isinstance(callable_t, ResolvedType)
                return callable_t.name
            return QualifiedName(result[1:])

        return result

    def report_error(self, error: ParserError):
        if isinstance(error, NameResolutionError):
            if error.name[0] in ["PyCapsule"]:
                return
        super().report_error(error)


class FixRedundantBuiltinsAnnotation(IParser):
    def handle_attribute(self, path: QualifiedName, attr: Any) -> Attribute | None:
        result = super().handle_attribute(path, attr)
        if result is None:
            return None
        if attr is None or inspect.ismodule(attr):
            result.annotation = None
        return result


class FixMissingNoneHashFieldAnnotation(IParser):
    def handle_field(self, path: QualifiedName, field: Any) -> Field | None:
        result = super().handle_field(path, field)
        if result is None:
            return None
        if field is None and path[-1] == "__hash__":
            result.attribute.annotation = self.parse_annotation_str("typing.ClassVar[None]")
        return result


class FixPEP585CollectionNames(IParser):
    __typing_collection_names: set[Identifier] = set(
        Identifier(name)
        for name in (
            "Dict",
            "List",
            "Set",
            "Tuple",
            "FrozenSet",
            "Type",
        )
    )

    def parse_annotation_str(self, annotation_str: str) -> ResolvedType | InvalidExpression | Value:
        result = super().parse_annotation_str(annotation_str)
        if not isinstance(result, ResolvedType) or len(result.name) != 2 or result.name[0] != "typing":
            return result

        word = result.name[1]
        if word in self.__typing_collection_names:
            result.name = QualifiedName.from_str(f"{word.lower()}")

        return result


class FixTypingTypeNames(IParser):
    __typing_names: set[Identifier] = set(
        Identifier(name)
        for name in (
            "Annotated",
            "Any",
            "Buffer",
            "Callable",
            "Dict",
            "ItemsView",
            "Iterable",
            "Iterator",
            "KeysView",
            "List",
            "Literal",
            "Optional",
            "Sequence",
            "Set",
            "Tuple",
            "Union",
            "ValuesView",
            "buffer",
            "iterable",
            "iterator",
            "sequence",
        )
    )
    __typing_extensions_names: set[Identifier] = set(
        Identifier(name)
        for name in (
            "buffer",
            "Buffer",
        )
    )

    def __init__(self):
        super().__init__()
        if sys.version_info < (3, 9):
            self.__typing_extensions_names.add(Identifier("Annotated"))
        if sys.version_info < (3, 8):
            self.__typing_extensions_names.add(Identifier("Literal"))

    def parse_annotation_str(self, annotation_str: str) -> ResolvedType | InvalidExpression | Value:
        result = super().parse_annotation_str(annotation_str)
        return self._parse_annotation_str(result)

    def _parse_annotation_str(
        self, result: ResolvedType | InvalidExpression | Value
    ) -> ResolvedType | InvalidExpression | Value:
        if not isinstance(result, ResolvedType):
            return result

        result.parameters = (
            [self._parse_annotation_str(p) for p in result.parameters] if result.parameters is not None else None
        )

        if len(result.name) != 1:
            return result

        word = result.name[0]
        if word in self.__typing_names:
            package = "typing"
            if word in self.__typing_extensions_names:
                package = "typing_extensions"
            result.name = QualifiedName.from_str(f"{package}.{word[0].upper()}{word[1:]}")
        if word == "function" and result.parameters is None:
            result.name = QualifiedName.from_str("typing.Callable")
        if word in ("object", "handle") and result.parameters is None:
            result.name = QualifiedName.from_str("typing.Any")

        return result


class FixCurrentModulePrefixInTypeNames(IParser):
    def __init__(self):
        super().__init__()
        self.__current_module: QualifiedName = QualifiedName()

    def handle_alias(self, path: QualifiedName, origin: Any) -> Alias | None:
        result = super().handle_alias(path, origin)
        if result is None:
            return None
        result.origin = self._strip_current_module(result.origin)
        return result

    def handle_attribute(self, path: QualifiedName, attr: Any) -> Attribute | None:
        result = super().handle_attribute(path, attr)
        if result is None:
            return None
        if isinstance(result.annotation, ResolvedType):
            result.annotation.name = self._strip_current_module(result.annotation.name)
        return result

    def handle_module(self, path: QualifiedName, module: types.ModuleType) -> Module | None:
        tmp = self.__current_module
        self.__current_module = path
        result = super().handle_module(path, module)
        self.__current_module = tmp
        return result

    def handle_type(self, type_: type) -> QualifiedName:
        result = super().handle_type(type_)
        return self._strip_current_module(result)

    def handle_value(self, value: Any) -> Value:
        result = super().handle_value(value)
        if inspect.isroutine(value):
            result.repr = str(self._strip_current_module(QualifiedName.from_str(result.repr)))
        return result

    def parse_annotation_str(self, annotation_str: str) -> ResolvedType | InvalidExpression | Value:
        result = super().parse_annotation_str(annotation_str)
        if isinstance(result, ResolvedType):
            result.name = self._strip_current_module(result.name)
        return result

    def _strip_current_module(self, name: QualifiedName) -> QualifiedName:
        if name[: len(self.__current_module)] == self.__current_module:
            return QualifiedName(name[len(self.__current_module) :])
        return name


class FixValueReprRandomAddress(IParser):
    _pattern = re.compile(
        r"<(?P<name>[\w.]+) object "
        r"(?P<capsule>\w+\s)*at "
        r"(?P<address>0x[a-fA-F0-9]+)>"
    )

    def handle_value(self, value: Any) -> Value:
        result = super().handle_value(value)
        result.repr = self._pattern.sub(r"<\g<name> object>", result.repr)
        return result


class FixNumpyArrayDimAnnotation(IParser):
    __array_names: set[QualifiedName] = {
        QualifiedName.from_str("numpy.ndarray"),
        *(
            QualifiedName.from_str(f"scipy.sparse.{storage}_{arr}")
            for storage in ["bsr", "coo", "csr", "csc", "dia", "dok", "lil"]
            for arr in ["array", "matrix"]
        ),
    }
    __annotated_name = QualifiedName.from_str("Annotated")
    numpy_primitive_types: set[QualifiedName] = set(
        map(
            QualifiedName.from_str,
            (
                "bool",
                *map(
                    lambda name: f"numpy.{name}",
                    (
                        "uint8",
                        "int8",
                        "uint16",
                        "int16",
                        "uint32",
                        "int32",
                        "uint64",
                        "int64",
                        "float16",
                        "float32",
                        "float64",
                        "complex32",
                        "complex64",
                        "longcomplex",
                    ),
                ),
            ),
        )
    )

    __DIM_VARS = ["n", "m"]

    def parse_annotation_str(self, annotation_str: str) -> ResolvedType | InvalidExpression | Value:
        result = super().parse_annotation_str(annotation_str)
        if (
            not isinstance(result, ResolvedType)
            or result.name not in self.__array_names
            or result.parameters is None
            or len(result.parameters) == 0
        ):
            return result

        scalar_with_dims = result.parameters[0]
        flags = result.parameters[1:]

        if not isinstance(scalar_with_dims, ResolvedType) or scalar_with_dims.name not in self.numpy_primitive_types:
            return result

        result = ResolvedType(
            name=self.__annotated_name,
            parameters=[
                self.parse_annotation_str(str(result.name)),
                ResolvedType(scalar_with_dims.name),
            ],
        )

        assert result.parameters is not None
        if scalar_with_dims.parameters is not None and len(scalar_with_dims.parameters) >= 0:
            dims = self.__to_dims(scalar_with_dims.parameters)
            if dims is not None and len(dims) > 0:
                result.parameters += [self.handle_value(self.__wrap_with_size_helper(dims))]

        result.parameters += flags

        return result

    def __wrap_with_size_helper(self, dims: list[int | str]) -> FixedSize | DynamicSize:
        if all(isinstance(d, int) for d in dims):
            return_t = FixedSize
        else:
            return_t = DynamicSize

        self.handle_type(return_t)
        return return_t(*dims)

    def __to_dims(self, dimensions: Sequence[ResolvedType | Value | InvalidExpression]) -> list[int | str] | None:
        result = []
        for dim_param in dimensions:
            if isinstance(dim_param, Value):
                try:
                    dim = int(dim_param.repr)
                except ValueError:
                    return None
            elif isinstance(dim_param, ResolvedType):
                dim = str(dim_param)
                if dim not in self.__DIM_VARS:
                    return None
            else:
                return None
            result.append(dim)
        return result

    def report_error(self, error: ParserError) -> None:
        if (
            isinstance(error, NameResolutionError)
            and len(error.name) == 1
            and len(error.name[0]) == 1
            and error.name[0] in self.__DIM_VARS
        ):
            return
        super().report_error(error)


class FixNumpyArrayDimTypeVar(IParser):
    __array_names: set[QualifiedName] = {QualifiedName.from_str("numpy.ndarray")}
    numpy_primitive_types = FixNumpyArrayDimAnnotation.numpy_primitive_types

    __DIM_VARS: set[str] = set()

    def handle_module(self, path: QualifiedName, module: types.ModuleType) -> Module | None:
        result = super().handle_module(path, module)
        if result is None:
            return None

        if self.__DIM_VARS:
            result.imports.add(Import(name=None, origin=QualifiedName.from_str("typing")))

            for name in self.__DIM_VARS:
                result.type_vars.append(
                    TypeVar_(
                        name=Identifier(name),
                        bound=self.parse_annotation_str("int"),
                    ),
                )

        self.__DIM_VARS.clear()

        return result

    def parse_annotation_str(self, annotation_str: str) -> ResolvedType | InvalidExpression | Value:
        result = super().parse_annotation_str(annotation_str)

        if not isinstance(result, ResolvedType):
            return result

        if len(result.name) == 1 and len(result.name[0]) == 1:
            result.name = QualifiedName.from_str(result.name[0].upper())
            self.__DIM_VARS.add(result.name[0])

        if result.name not in self.__array_names:
            return result

        if result.parameters is None or len(result.parameters) == 0:
            result.parameters = [
                self.parse_annotation_str("Any"),
                ResolvedType(
                    name=QualifiedName.from_str("numpy.dtype"),
                    parameters=[self.parse_annotation_str("Any")],
                ),
            ]
            return result

        scalar_with_dims = result.parameters[0]

        if not isinstance(scalar_with_dims, ResolvedType) or scalar_with_dims.name not in self.numpy_primitive_types:
            return result

        name = scalar_with_dims.name
        if str(name) == "bool":
            name = QualifiedName.from_str("numpy.bool_")
        dtype = ResolvedType(
            name=QualifiedName.from_str("numpy.dtype"),
            parameters=[ResolvedType(name=name)],
        )

        shape = self.parse_annotation_str("Any")
        if scalar_with_dims.parameters is not None and len(scalar_with_dims.parameters) > 0:
            dims = self.__to_dims(scalar_with_dims.parameters)
            if dims is not None:
                shape = self.parse_annotation_str("Tuple")
                assert isinstance(shape, ResolvedType)
                shape.parameters = []
                for dim in dims:
                    if isinstance(dim, int):
                        literal_dim = self.parse_annotation_str("Literal")
                        assert isinstance(literal_dim, ResolvedType)
                        literal_dim.parameters = [Value(repr=str(dim))]
                        shape.parameters.append(literal_dim)
                    else:
                        shape.parameters.append(ResolvedType(name=QualifiedName.from_str(dim)))

        result.parameters = [shape, dtype]
        return result

    def __to_dims(self, dimensions: Sequence[ResolvedType | Value | InvalidExpression]) -> list[int | str] | None:
        result: list[int | str] = []
        for dim_param in dimensions:
            if isinstance(dim_param, Value):
                try:
                    dim = int(dim_param.repr)
                except ValueError:
                    return None
            elif isinstance(dim_param, ResolvedType):
                dim = str(dim_param)
            else:
                return None
            result.append(dim)
        return result

    def report_error(self, error: ParserError) -> None:
        if isinstance(error, NameResolutionError) and len(error.name) == 1 and error.name[0] in self.__DIM_VARS:
            return
        super().report_error(error)


class FixNumpyArrayRemoveParameters(IParser):
    __ndarray_name = QualifiedName.from_str("numpy.ndarray")

    def parse_annotation_str(self, annotation_str: str) -> ResolvedType | InvalidExpression | Value:
        result = super().parse_annotation_str(annotation_str)
        if isinstance(result, ResolvedType) and result.name == self.__ndarray_name:
            result.parameters = None
        return result


class FixScipyTypeArguments(IParser):
    def parse_annotation_str(self, annotation_str: str) -> ResolvedType | InvalidExpression | Value:
        result = super().parse_annotation_str(annotation_str)

        if not isinstance(result, ResolvedType):
            return result

        if result.name[:2] == ("scipy", "sparse"):
            result.parameters = None

        return result


class FixNumpyDtype(IParser):
    __numpy_dtype = QualifiedName.from_str("numpy.dtype")

    def parse_annotation_str(self, annotation_str: str) -> ResolvedType | InvalidExpression | Value:
        result = super().parse_annotation_str(annotation_str)

        if not isinstance(result, ResolvedType) or result.parameters:
            return result

        if result.name[:1] == ("dtype",) or result.name[:2] == ("numpy", "dtype"):
            result.name = self.__numpy_dtype
            result.parameters = [self.parse_annotation_str("Any")]

        return result


class FixNumpyArrayFlags(IParser):
    __ndarray_name = QualifiedName.from_str("numpy.ndarray")
    __flags: set[QualifiedName] = {
        QualifiedName.from_str("flags.writeable"),
        QualifiedName.from_str("flags.c_contiguous"),
        QualifiedName.from_str("flags.f_contiguous"),
    }

    def parse_annotation_str(self, annotation_str: str) -> ResolvedType | InvalidExpression | Value:
        result = super().parse_annotation_str(annotation_str)
        if isinstance(result, ResolvedType) and result.name == self.__ndarray_name:
            if result.parameters is not None:
                for param in result.parameters:
                    if isinstance(param, ResolvedType) and param.name in self.__flags:
                        param.name = QualifiedName.from_str(f"numpy.ndarray.{param.name}")

        return result

    def report_error(self, error: ParserError) -> None:
        if isinstance(error, NameResolutionError) and error.name in self.__flags:
            return
        super().report_error(error)


class FixRedundantMethodsFromBuiltinObject(IParser):
    def handle_method(self, path: QualifiedName, method: Any) -> list[Method]:
        result = super().handle_method(path, method)
        return [m for m in result if not (m.function.name == "__init__" and m.function.doc == object.__init__.__doc__)]


class ReplaceReadWritePropertyWithField(IParser):
    def handle_class_member(
        self, path: QualifiedName, class_: type, obj: Any
    ) -> Docstring | Alias | Class | list[Method] | Field | Property | None:
        result = super().handle_class_member(path, class_, obj)
        if isinstance(result, Property):
            if (
                result.doc is None
                and result.getter is not None
                and result.setter is not None
                and len(result.getter.args) == 1
                and len(result.setter.args) == 2
                and result.getter.doc is None
                and result.setter.doc is None
                and result.getter.returns == result.setter.args[1].annotation
            ):
                return Field(
                    attribute=Attribute(name=result.name, annotation=result.getter.returns, value=None),
                    modifier=None,
                )
        return result


class FixMissingFixedSizeImport(IParser):
    def parse_annotation_str(self, annotation_str: str) -> ResolvedType | InvalidExpression | Value:
        result = super().parse_annotation_str(annotation_str)
        if isinstance(result, Value) and result.repr.startswith("FixedSize(") and result.repr.endswith(")"):
            try:
                dimensions = map(
                    int,
                    result.repr[len("FixedSize(") : -len(")")].split(","),
                )
            except ValueError:
                pass
            else:
                self.handle_type(FixedSize)
                return self.handle_value(FixedSize(*dimensions))
        return result


class FixMissingEnumMembersAnnotation(IParser):
    __class_var_dict = ResolvedType(
        name=QualifiedName.from_str("typing.ClassVar"),
        parameters=[ResolvedType(name=QualifiedName.from_str("dict"))],
    )

    def handle_field(self, path: QualifiedName, field: Any) -> Field | None:
        result = super().handle_field(path, field)
        if result is None:
            return None
        if (
            path[-1] == "__members__"
            and isinstance(field, dict)
            and result.attribute.annotation == self.__class_var_dict
        ):
            assert isinstance(result.attribute.annotation, ResolvedType)
            dict_type = self._guess_dict_type(field)
            if dict_type is not None:
                result.attribute.annotation.parameters = [dict_type]
        return result

    def _guess_dict_type(self, d: dict) -> ResolvedType | None:
        if len(d) == 0:
            return None
        key_types = set()
        value_types = set()
        for key, value in d.items():
            key_types.add(self.handle_type(type(key)))
            value_types.add(self.handle_type(type(value)))
        if len(key_types) == 1:
            key_type = [ResolvedType(name=t) for t in key_types][0]
        else:
            union_t = self.parse_annotation_str("typing.Union")
            assert isinstance(union_t, ResolvedType)
            key_type = ResolvedType(name=union_t.name, parameters=[ResolvedType(name=t) for t in key_types])
        if len(value_types) == 1:
            value_type = [ResolvedType(name=t) for t in value_types][0]
        else:
            union_t = self.parse_annotation_str("typing.Union")
            assert isinstance(union_t, ResolvedType)
            value_type = ResolvedType(
                name=union_t.name,
                parameters=[ResolvedType(name=t) for t in value_types],
            )
        dict_t = self.parse_annotation_str("typing.Dict")
        assert isinstance(dict_t, ResolvedType)
        return ResolvedType(
            name=dict_t.name,
            parameters=[key_type, value_type],
        )


class FixPybind11EnumStrDoc(IParser):
    def handle_class_member(
        self, path: QualifiedName, class_: type, obj: Any
    ) -> Docstring | Alias | Class | list[Method] | Field | Property | None:
        result = super().handle_class_member(path, class_, obj)
        if not isinstance(result, list) or not hasattr(class_, "__members__"):
            return result
        for method in result:
            assert isinstance(method, Method)
            if method.function.name != "__str__" or method.function.doc != "name(self: handle) -> str\n":
                continue
            method.function.args = [
                Argument(
                    name=Identifier("self"),
                )
            ]
            method.function.returns = ResolvedType(name=QualifiedName.from_str("str"))
            method.modifier = None
            method.function.doc = None
        return result


class OverridePrintSafeValues(IParser):
    _print_safe_values: re.Pattern | None

    def __init__(self):
        super().__init__()
        self._print_safe_values = None

    def set_print_safe_value_pattern(self, pattern: re.Pattern):
        self._print_safe_values = pattern

    def parse_value_str(self, value: str) -> Value | InvalidExpression:
        result = super().parse_value_str(value)
        if (
            self._print_safe_values is not None
            and isinstance(result, Value)
            and not result.is_print_safe
            and self._print_safe_values.match(result.repr) is not None
        ):
            result.is_print_safe = True
        return result


class RewritePybind11EnumValueRepr(IParser):
    _pybind11_enum_pattern = re.compile(r"<(?P<enum>\w+(\.\w+)+): (?P<value>-?\d+)>")
    _unknown_enum_classes: set[str] = set()

    def __init__(self):
        super().__init__()
        self._pybind11_enum_locations: dict[re.Pattern, str] = {}

    def set_pybind11_enum_locations(self, locations: dict[re.Pattern, str]):
        self._pybind11_enum_locations = locations

    def parse_value_str(self, value: str) -> Value | InvalidExpression:
        value = value.strip()
        match = self._pybind11_enum_pattern.match(value)
        if match is not None:
            enum_qual_name = match.group("enum")
            enum_class_str, entry = enum_qual_name.rsplit(".", maxsplit=1)
            for pattern, prefix in self._pybind11_enum_locations.items():
                if pattern.match(enum_class_str) is None:
                    continue
                enum_class = self.parse_annotation_str(f"{prefix}.{enum_class_str}")
                if isinstance(enum_class, ResolvedType):
                    return Value(repr=f"{enum_class.name}.{entry}", is_print_safe=True)
        return super().parse_value_str(value)

    def report_error(self, error: ParserError) -> None:
        if isinstance(error, InvalidExpressionError):
            match = self._pybind11_enum_pattern.match(error.expression)
            if match is not None:
                enum_qual_name = match.group("enum")
                enum_class_str, entry = enum_qual_name.rsplit(".", maxsplit=1)
                self._unknown_enum_classes.add(enum_class_str)
        super().report_error(error)

    def finalize(self):
        if self._unknown_enum_classes:
            logger.warning(
                "Enum-like str representations were found with no "
                "matching mapping to the enum class location.\n"
                "Use `--enum-class-locations` to specify "
                "full path to the following enum(s):\n" + "\n".join(f" - {c}" for c in self._unknown_enum_classes)
            )
        super().finalize()


def indent_lines(lines: list[str], by=4) -> list[str]:
    return [" " * by + line for line in lines]


class Printer:
    def __init__(self, invalid_expr_as_ellipses: bool):
        self.invalid_expr_as_ellipses = invalid_expr_as_ellipses

    def print_alias(self, alias: Alias) -> list[str]:
        return [f"{alias.name} = {alias.origin}"]

    def print_attribute(self, attr: Attribute) -> list[str]:
        parts = [
            f"{attr.name}",
        ]
        if attr.annotation is not None:
            parts.append(f": {self.print_annotation(attr.annotation)}")

        if attr.value is not None and attr.value.is_print_safe:
            parts.append(f" = {self.print_value(attr.value)}")
        else:
            if attr.annotation is None:
                parts.append(" = ...")
            if attr.value is not None:
                parts.append(f"  # value = {self.print_value(attr.value)}")

        return ["".join(parts)]

    def print_argument(self, arg: Argument) -> str:
        parts = []
        if arg.variadic:
            parts += ["*"]
        if arg.kw_variadic:
            parts += ["**"]
        parts.append(f"{arg.name}")
        if arg.annotation is not None:
            parts.append(f": {self.print_annotation(arg.annotation)}")
        if isinstance(arg.default, Value):
            if arg.default.is_print_safe:
                parts.append(f" = {self.print_value(arg.default)}")
            else:
                parts.append(" = ...")
        elif isinstance(arg.default, InvalidExpression):
            parts.append(f" = {self.print_invalid_exp(arg.default)}")

        return "".join(parts)

    def print_class(self, class_: Class) -> list[str]:
        if class_.bases:
            base_list = "(" + ", ".join(str(base) for base in class_.bases) + ")"
        else:
            base_list = ""
        return [
            f"class {class_.name}{base_list}:",
            *indent_lines(self.print_class_body(class_)),
        ]

    def print_type_var(self, type_var: TypeVar_) -> list[str]:
        return [str(type_var)]

    def print_class_body(self, class_: Class) -> list[str]:
        result = []
        if class_.doc is not None:
            result.extend(self.print_docstring(class_.doc))

        for sub_class in sorted(class_.classes, key=lambda c: c.name):
            result.extend(self.print_class(sub_class))

        modifier_order: dict[Modifier, int] = {
            "static": 0,
            "class": 1,
            None: 2,
        }
        for field in sorted(class_.fields, key=lambda f: (modifier_order[f.modifier], f.attribute.name)):
            result.extend(self.print_field(field))

        for alias in sorted(class_.aliases, key=lambda a: a.name):
            result.extend(self.print_alias(alias))

        for method in sorted(class_.methods, key=lambda m: (modifier_order[m.modifier], m.function.name)):
            result.extend(self.print_method(method))

        for prop in sorted(class_.properties, key=lambda p: p.name):
            result.extend(self.print_property(prop))

        if not result:
            result = ["pass"]

        return result

    def print_docstring(self, doc: Docstring) -> list[str]:
        return [
            '"""',
            *(line.replace("\\", r"\\").replace('"""', r"\"\"\"") for line in doc.splitlines()),
            '"""',
        ]

    def print_field(self, field: Field) -> list[str]:
        return self.print_attribute(field.attribute)

    def print_function(self, func: Function) -> list[str]:
        pos_only = False
        kw_only = False

        args = []
        for arg in func.args:
            if arg.variadic:
                pos_only = True
                kw_only = True
            if not pos_only and not arg.pos_only:
                pos_only = True
                if sys.version_info >= (3, 8):
                    args.append("/")
            if not kw_only and arg.kw_only:
                kw_only = True
                args.append("*")
            args.append(self.print_argument(arg))
        if len(args) > 0 and args[0] == "/":
            args = args[1:]
        signature = [
            f"def {func.name}(",
            ", ".join(args),
            ")",
        ]

        if func.returns is not None:
            signature.append(f" -> {self.print_annotation(func.returns)}")
        signature.append(":")

        result: list[str] = [
            *(f"@{decorator}" for decorator in func.decorators),
            "".join(signature),
        ]

        if func.doc is not None:
            body = self.print_docstring(func.doc)
        else:
            body = ["..."]

        result.extend(indent_lines(body))

        return result

    def print_submodule_import(self, name: Identifier) -> list[str]:
        return [f"from . import {name}"]

    def print_import(self, import_: Import) -> list[str]:
        parent = str(import_.origin.parent)
        if import_.name is None:
            return [f"import {import_.origin}"]

        if len(parent) == 0:
            return [f"import {import_.origin} as {import_.name}"]

        result = f"from {parent} import {import_.origin[-1]}"
        if import_.name != import_.origin[-1]:
            result += f" as {import_.name}"
        return [result]

    def print_method(self, method: Method) -> list[str]:
        result = []
        if method.modifier == "static":
            result += ["@staticmethod"]
        elif method.modifier == "class":
            result += ["@classmethod"]
        elif method.modifier is None:
            pass
        else:
            raise RuntimeError()
        result.extend(self.print_function(method.function))
        return result

    def print_module(self, module: Module) -> list[str]:
        result = []

        if module.doc is not None:
            result.extend(self.print_docstring(module.doc))

        for import_ in sorted(module.imports, key=lambda x: x.origin):
            result.extend(self.print_import(import_))

        for sub_module in module.sub_modules:
            result.extend(self.print_submodule_import(sub_module.name))

        for attr in sorted(module.attributes, key=lambda a: a.name):
            if attr.name == "__all__":
                result.extend(self.print_attribute(attr))
                break

        for type_var in sorted(module.type_vars, key=lambda t: t.name):
            result.extend(self.print_type_var(type_var))

        for class_ in sorted(module.classes, key=lambda c: c.name):
            result.extend(self.print_class(class_))

        for func in sorted(module.functions, key=lambda f: f.name):
            result.extend(self.print_function(func))

        for attr in sorted(module.attributes, key=lambda a: a.name):
            if attr.name != "__all__":
                result.extend(self.print_attribute(attr))

        for alias in module.aliases:
            result.extend(self.print_alias(alias))

        return result

    def print_property(self, prop: Property) -> list[str]:
        if not prop.getter:
            return []

        result = []

        result.extend(
            [
                "@property",
                *self.print_function(
                    dataclasses.replace(
                        prop.getter,
                        name=prop.name,
                        doc=prop.doc if prop.doc is not None else prop.getter.doc,
                    )
                ),
            ]
        )
        if prop.setter:
            result.extend(
                [
                    f"@{prop.name}.setter",
                    *self.print_function(
                        dataclasses.replace(
                            prop.setter,
                            name=prop.name,
                            doc=None if prop.doc is not None else prop.setter.doc,
                        )
                    ),
                ]
            )

        return result

    def print_value(self, value: Value) -> str:
        split = value.repr.split("\n", 1)
        if len(split) == 1:
            return split[0]
        else:
            return split[0] + "..."

    def print_type(self, type_: ResolvedType) -> str:
        if str(type_.name) == "typing.Optional" and type_.parameters is not None and len(type_.parameters) == 1:
            return f"{self.print_annotation(type_.parameters[0])} | None"
        if str(type_.name) == "typing.Union" and type_.parameters is not None:
            return " | ".join(self.print_annotation(p) for p in type_.parameters)
        if type_.parameters:
            param_str = "[" + ", ".join(self.print_annotation(p) for p in type_.parameters) + "]"
        else:
            param_str = ""
        return f"{type_.name}{param_str}"

    def print_annotation(self, annotation: Annotation) -> str:
        if isinstance(annotation, ResolvedType):
            return self.print_type(annotation)
        elif isinstance(annotation, Value):
            return self.print_value(annotation)
        elif isinstance(annotation, InvalidExpression):
            return self.print_invalid_exp(annotation)
        else:
            raise AssertionError()

    def print_invalid_exp(self, invalid_expr: InvalidExpression) -> str:
        if self.invalid_expr_as_ellipses:
            return "..."
        return invalid_expr.text


class Writer:
    def __init__(self, stub_ext: str = "pyi"):
        self.stub_ext: str = stub_ext

    def write_module(self, module: Module, printer: Printer, to: Path, sub_dir: Path | None = None):
        assert to.exists()
        assert to.is_dir()
        if module.sub_modules or module.is_package or sub_dir is not None:
            if sub_dir is None:
                sub_dir = Path(module.name)
            module_dir = to / sub_dir
            module_dir.mkdir(exist_ok=True)
            module_file = module_dir / f"__init__.{self.stub_ext}"
        else:
            module_file = to / f"{module.name}.{self.stub_ext}"
        with open(module_file, "w", encoding="utf-8") as f:
            f.writelines(line + "\n" for line in printer.print_module(module))

        for sub_module in module.sub_modules:
            self.write_module(sub_module, printer, to=module_dir)


class CLIArgs(Namespace):
    output_dir: str
    root_suffix: str
    ignore_invalid_expressions: re.Pattern | None
    ignore_invalid_identifiers: re.Pattern | None
    ignore_unresolved_names: re.Pattern | None
    ignore_all_errors: bool
    enum_class_locations: list[tuple[re.Pattern, str]]
    numpy_array_wrap_with_annotated: bool
    numpy_array_use_type_var: bool
    numpy_array_remove_parameters: bool
    print_invalid_expressions_as_is: bool
    print_safe_value_reprs: re.Pattern | None
    exit_code: bool
    dry_run: bool
    stub_extension: str
    module_name: str


def arg_parser() -> ArgumentParser:
    def regex(pattern_str: str) -> re.Pattern:
        try:
            return re.compile(pattern_str)
        except re.error as e:
            raise ValueError(f"Invalid REGEX pattern: {e}")

    def regex_colon_path(regex_path: str) -> tuple[re.Pattern, str]:
        pattern_str, path = regex_path.rsplit(":", maxsplit=1)
        if any(not part.isidentifier() for part in path.split(".")):
            raise ValueError(f"Invalid PATH: {path}")
        return regex(pattern_str), path

    parser = ArgumentParser(prog="pybind11-stubgen", description="Generates stubs for specified modules")
    parser.add_argument(
        "-o",
        "--output-dir",
        help="The root directory for output stubs",
        default="./stubs",
    )
    parser.add_argument(
        "--root-suffix",
        type=str,
        default=None,
        dest="root_suffix",
        help="Top-level module directory suffix",
    )

    parser.add_argument(
        "--ignore-invalid-expressions",
        metavar="REGEX",
        default=None,
        type=regex,
        help="Ignore invalid expressions matching REGEX",
    )
    parser.add_argument(
        "--ignore-invalid-identifiers",
        metavar="REGEX",
        default=None,
        type=regex,
        help="Ignore invalid identifiers matching REGEX",
    )

    parser.add_argument(
        "--ignore-unresolved-names",
        metavar="REGEX",
        default=None,
        type=regex,
        help="Ignore unresolved names matching REGEX",
    )

    parser.add_argument(
        "--ignore-all-errors",
        default=False,
        action="store_true",
        help="Ignore all errors during module parsing",
    )

    parser.add_argument(
        "--enum-class-locations",
        dest="enum_class_locations",
        metavar="REGEX:LOC",
        action="append",
        default=[],
        type=regex_colon_path,
        help="Locations of enum classes in "
        "<enum-class-name-regex>:<path-to-class> format. "
        "Example: `MyEnum:foo.bar.Baz`",
    )

    numpy_array_fix = parser.add_mutually_exclusive_group()
    numpy_array_fix.add_argument(
        "--numpy-array-wrap-with-annotated",
        default=False,
        action="store_true",
        help="Replace numpy/scipy arrays of "
        "'ARRAY_T[TYPE, [*DIMS], *FLAGS]' format with "
        "'Annotated[ARRAY_T, TYPE, FixedSize|DynamicSize(*DIMS), *FLAGS]'",
    )
    numpy_array_fix.add_argument(
        "--numpy-array-use-type-var",
        default=False,
        action="store_true",
        help="Replace 'numpy.ndarray[numpy.float32[m, 1]]' with "
        "'numpy.ndarray[tuple[M, typing.Literal[1]], numpy.dtype[numpy.float32]]'",
    )

    numpy_array_fix.add_argument(
        "--numpy-array-remove-parameters",
        default=False,
        action="store_true",
        help="Replace 'numpy.ndarray[...]' with 'numpy.ndarray'",
    )

    parser.add_argument(
        "--print-invalid-expressions-as-is",
        default=False,
        action="store_true",
        help="Suppress the replacement with '...' of invalid expressionsfound in annotations",
    )

    parser.add_argument(
        "--print-safe-value-reprs",
        metavar="REGEX",
        default=None,
        type=regex,
        help="Override the print-safe check for values matching REGEX",
    )

    parser.add_argument(
        "--exit-code",
        action="store_true",
        dest="exit_code",
        help="On error exits with 1 and skips stub generation",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        dest="dry_run",
        help="Don't write stubs. Parses module and report errors",
    )

    parser.add_argument(
        "--stub-extension",
        type=str,
        default="pyi",
        metavar="EXT",
        choices=["pyi", "py"],
        help="The file extension of the generated stubs. Must be 'pyi' (default) or 'py'",
    )

    parser.add_argument(
        "module_name",
        metavar="MODULE_NAME",
        type=str,
        help="module name",
    )

    return parser


def stub_parser_from_args(args: CLIArgs) -> IParser:
    error_handlers_top: list[type] = [
        LoggerData,
        *([IgnoreAllErrors] if args.ignore_all_errors else []),
        *([IgnoreInvalidIdentifierErrors] if args.ignore_invalid_identifiers else []),
        *([IgnoreInvalidExpressionErrors] if args.ignore_invalid_expressions else []),
        *([IgnoreUnresolvedNameErrors] if args.ignore_unresolved_names else []),
    ]
    error_handlers_bottom: list[type] = [
        LogErrors,
        SuggestCxxSignatureFix,
        *([TerminateOnFatalErrors] if args.exit_code else []),
    ]

    numpy_fixes: list[type] = [
        *([FixNumpyArrayDimAnnotation] if args.numpy_array_wrap_with_annotated else []),
        *([FixNumpyArrayDimTypeVar] if args.numpy_array_use_type_var else []),
        *([FixNumpyArrayRemoveParameters] if args.numpy_array_remove_parameters else []),
    ]

    # Parser composition follows after mixin class definitions

    class ParserDispatchMixin(IParser):
        def handle_class(self, path: QualifiedName, class_: type) -> Class | None:
            base_classes = class_.__bases__
            result = Class(name=path[-1], bases=self.handle_bases(path, base_classes))
            for name, member in inspect.getmembers(class_):
                obj = self.handle_class_member(QualifiedName([*path, Identifier(name)]), class_, member)
                if isinstance(obj, Docstring):
                    result.doc = obj
                elif isinstance(obj, Alias):
                    result.aliases.append(obj)
                elif isinstance(obj, Field):
                    result.fields.append(obj)
                elif isinstance(obj, Class):
                    result.classes.append(obj)
                elif isinstance(obj, list):
                    result.methods.extend(obj)
                elif isinstance(obj, Property):
                    result.properties.append(obj)
                elif obj is None:
                    pass
                else:
                    raise AssertionError()
            return result

        def handle_class_member(
            self, path: QualifiedName, class_: type, member: Any
        ) -> Docstring | Alias | Class | list[Method] | Field | Property | None:
            if inspect.isroutine(member):
                return self.handle_method(path, member)
            if self._is_alias(path, member):
                return self.handle_alias(path, member)
            if inspect.isclass(member):
                return self.handle_class(path, member)
            if self._is_descriptor(member):
                return self.handle_property(path, member)
            if path[-1] == "__doc__":
                return self.handle_docstring(path, member)
            return self.handle_field(path, member)

        def handle_module(self, path: QualifiedName, module: types.ModuleType) -> Module | None:
            result = Module(name=path[-1])
            for name, member in inspect.getmembers(module):
                obj = self.handle_module_member(QualifiedName([*path, Identifier(name)]), module, member)
                if isinstance(obj, Docstring):
                    result.doc = obj
                elif isinstance(obj, Import):
                    result.imports.add(obj)
                elif isinstance(obj, Alias):
                    result.aliases.append(obj)
                elif isinstance(obj, Class):
                    result.classes.append(obj)
                elif isinstance(obj, list):
                    result.functions.extend(obj)
                elif isinstance(obj, Module):
                    result.sub_modules.append(obj)
                elif isinstance(obj, Attribute):
                    result.attributes.append(obj)
                elif isinstance(obj, TypeVar_):
                    result.type_vars.append(obj)
                elif obj is None:
                    if name == "__path__":
                        result.is_package = True
                else:
                    raise AssertionError()

            return result

        def handle_module_member(
            self, path: QualifiedName, module: types.ModuleType, member: Any
        ) -> Docstring | Import | Alias | Class | list[Function] | Attribute | Module | None:
            member_module = self._get_value_parent_module_name(member)
            if member_module is not None and member_module != module.__name__ or path[-1] == "annotations":
                return self.handle_import(path, member)
            if self._is_alias(path, member):
                return self.handle_alias(path, member)
            if inspect.isroutine(member):
                return self.handle_function(path, member)
            if inspect.isclass(member):
                return self.handle_class(path, member)
            if inspect.ismodule(member):
                return self.handle_module(path, member)
            if path[-1] == "__doc__":
                return self.handle_docstring(path, member)
            return self.handle_attribute(path, member)

        def _get_value_parent_module_name(self, obj: Any) -> str | None:
            if inspect.ismodule(obj):
                return obj.__name__.rsplit(".", 1)[0]
            if inspect.isclass(obj) or inspect.isroutine(obj):
                return getattr(obj, "__module__", None)
            return None

        def _is_alias(self, path: QualifiedName, member: Any):
            if (inspect.isroutine(member) or inspect.isclass(member)) and path[-1] != member.__name__:
                return True
            if inspect.ismodule(member) and member.__name__ != str(path):
                return True
            return False

        def _is_descriptor(self, member: Any) -> bool:
            return hasattr(member, "__get__") or hasattr(member, "__set__")

        def report_error(self, error: ParserError) -> None:
            if isinstance(error, NameResolutionError):
                if error.name[0] == "pybind11_builtins":
                    return
            if isinstance(error, InvalidIdentifierError):
                name = error.name
                if (
                    name.startswith("ItemsView[")
                    and name.endswith("]")
                    or name.startswith("KeysView[")
                    and name.endswith("]")
                    or name.startswith("ValuesView[")
                    and name.endswith("]")
                ):
                    return

            super().report_error(error)

        def handle_bases(self, path: QualifiedName, bases: tuple[type, ...]) -> list[QualifiedName]:
            result = []
            for base in super().handle_bases(path, bases):
                if base[0] == "pybind11_builtins":
                    break
                result.append(base)
            return result

    class BaseParser(IParser):
        def handle_alias(self, path: QualifiedName, origin: Any) -> Alias | None:
            full_name = self._get_full_name(path, origin)
            if full_name is None:
                return None

            return Alias(
                name=path[-1],
                origin=full_name,
            )

        def handle_attribute(self, path: QualifiedName, value: Any) -> Attribute | None:
            return Attribute(
                name=path[-1],
                value=self.handle_value(value),
                annotation=ResolvedType(name=self.handle_type(type(value))),
            )

        def handle_bases(self, path: QualifiedName, bases: tuple[type, ...]) -> list[QualifiedName]:
            return [self.handle_type(type_) for type_ in bases if type_ is not object]

        def handle_docstring(self, path: QualifiedName, value: Any) -> Docstring | None:
            if isinstance(value, str):
                return Docstring(value)
            return None

        def handle_field(self, path: QualifiedName, value: Any) -> Field | None:
            attr = self.handle_attribute(path, value)
            if attr is None:
                return None
            result = Field(
                attribute=Attribute(
                    name=attr.name,
                    value=attr.value,
                ),
                modifier="static",
            )
            if attr.annotation is not None:
                class_var = self.parse_annotation_str("typing.ClassVar")
                assert isinstance(class_var, ResolvedType)
                result.attribute.annotation = ResolvedType(
                    name=class_var.name,
                    parameters=[attr.annotation],
                )
            return result

        def handle_function(self, path: QualifiedName, func: Any) -> list[Function]:
            doc: Docstring | None = Docstring(func.__doc__) if getattr(func, "__doc__", None) is not None else None
            func_name = Identifier(path[-1])

            try:
                (
                    args,
                    var_args,
                    var_kw,
                    defaults,
                    kw_only_args,
                    kw_only_defaults,
                    annotations,
                ) = inspect.getfullargspec(func)

                func_args: dict[str, Argument] = {arg_name: Argument(name=Identifier(arg_name)) for arg_name in args}
                func_args["return"] = Argument(
                    name=Identifier(),
                )
                if var_args is not None:
                    func_args[var_args] = Argument(name=Identifier(var_args), variadic=True)
                for arg_name in kw_only_args:
                    func_args[arg_name] = Argument(name=Identifier(arg_name), kw_only=True)
                if var_kw is not None:
                    func_args[var_kw] = Argument(name=Identifier(var_kw), kw_variadic=True)
                if defaults is not None:
                    for i, default in enumerate(defaults):
                        arg_name = args[len(args) - len(defaults) + i]
                        func_args[arg_name].default = self.handle_value(default)
                if kw_only_defaults is not None:
                    for arg_name, default in kw_only_defaults.items():
                        func_args[arg_name].default = self.handle_value(default)
                for arg_name, annotation in annotations.items():
                    if isinstance(annotation, str):
                        func_args[arg_name].annotation = self.parse_annotation_str(annotation)
                    elif not isinstance(annotation, type):
                        func_args[arg_name].annotation = self.handle_value(annotation)
                    elif self._is_generic_alias(annotation):
                        func_args[arg_name].annotation = self.parse_annotation_str(str(annotation))
                    else:
                        func_args[arg_name].annotation = ResolvedType(
                            name=self.handle_type(annotation),
                        )
                if "return" in func_args:
                    returns = func_args["return"].annotation
                else:
                    returns = None
                return [
                    Function(
                        name=func_name,
                        args=[arg for arg_name, arg in func_args.items() if arg_name != "return"],
                        returns=returns,
                        doc=doc,
                    )
                ]

            except TypeError:
                return [
                    Function(
                        name=func_name,
                        args=[
                            Argument(name=Identifier("args"), variadic=True),
                            Argument(name=Identifier("kwargs"), kw_variadic=True),
                        ],
                        doc=doc,
                    )
                ]

        def _is_generic_alias(self, annotation: type) -> bool:
            generic_alias_t: type | None = getattr(types, "GenericAlias", None)
            if generic_alias_t is None:
                return False
            return isinstance(annotation, generic_alias_t)

        def handle_import(self, path: QualifiedName, origin: Any) -> Import | None:
            full_name = self._get_full_name(path, origin)
            if full_name is None:
                return None

            return Import(
                path[-1],
                full_name,
            )

        def handle_method(self, path: QualifiedName, method: Any) -> list[Method]:
            return [
                Method(function=func, modifier=self._get_method_modifier(func.args))
                for func in self.handle_function(path, method)
            ]

        def handle_value(self, value: Any) -> Value:
            value_type = type(value)
            if value is None or value_type in (bool, int, str):
                return Value(repr=repr(value), is_print_safe=True)
            if value_type in (float, complex):
                try:
                    repr_str = repr(value)
                    eval(repr_str)
                    return Value(repr=repr_str, is_print_safe=True)
                except (SyntaxError, NameError):
                    pass
            if value_type in (list, tuple, set):
                assert isinstance(value, (list, tuple, set))
                if len(value) == 0:
                    return Value(repr=f"{value_type.__name__}()", is_print_safe=True)
                elements: list[Value] = [self.handle_value(el) for el in value]
                is_print_safe = all(el.is_print_safe for el in elements)
                left, right = {
                    list: "[]",
                    tuple: "()",
                    set: "{}",
                }[value_type]
                return Value(
                    repr="".join([left, ", ".join(el.repr for el in elements), right]),
                    is_print_safe=is_print_safe,
                )
            if value_type is dict:
                assert isinstance(value, dict)
                parts = []
                is_print_safe = True
                for k, v in value.items():
                    k_value = self.handle_value(k)
                    v_value = self.handle_value(v)
                    parts.append(f"{k_value.repr}: {v_value.repr}")
                    is_print_safe = is_print_safe and k_value.is_print_safe and v_value.is_print_safe

                return Value(
                    repr="".join(["{", ", ".join(parts), "}"]),
                    is_print_safe=is_print_safe,
                )
            if inspect.isroutine(value):
                module_name = getattr(value, "__module__", None)
                qual_name = getattr(value, "__qualname__", None)
                if (
                    module_name is not None
                    and "<" not in module_name
                    and qual_name is not None
                    and "<" not in qual_name
                ):
                    if module_name == "builtins":
                        repr_str = qual_name
                    else:
                        repr_str = f"{module_name}.{qual_name}"
                    return Value(repr=repr_str, is_print_safe=True)
            if inspect.isclass(value):
                return Value(repr=str(self.handle_type(value)), is_print_safe=True)
            if inspect.ismodule(value):
                return Value(repr=value.__name__, is_print_safe=True)
            return Value(repr=repr(value), is_print_safe=False)

        def handle_type(self, type_: type) -> QualifiedName:
            return QualifiedName(
                (
                    Identifier(part)
                    for part in (
                        *type_.__module__.split("."),
                        *type_.__qualname__.split("."),
                    )
                )
            )

        def parse_value_str(self, value: str) -> Value | InvalidExpression:
            return self._parse_expression_str(value)

        def report_error(self, error: ParserError):
            if isinstance(error, NameResolutionError):
                if error.name[0] == "module":
                    return
            super().report_error(error)

        def _get_method_modifier(self, args: list[Argument]) -> Modifier:
            if len(args) == 0:
                return "static"
            name = args[0].name
            if name == Identifier("self"):
                return None
            elif name == Identifier("cls"):
                return "class"
            else:
                return "static"

        def _get_full_name(self, path: QualifiedName, origin: Any) -> QualifiedName | None:
            if inspect.ismodule(origin):
                origin_full_name = origin.__name__
            else:
                module_name = getattr(origin, "__module__", None)
                if module_name == "__future__":
                    return None
                if module_name is None:
                    self.report_error(NameResolutionError(path))
                    return None
                qual_name = getattr(origin, "__qualname__", None)
                if qual_name is None:
                    self.report_error(NameResolutionError(path))
                    return None
                match = re.match(
                    r"(PyCapsule|pybind11_detail_function_record_[_a-zA-Z0-9]+)\.",
                    qual_name,
                )
                if match:
                    qual_name = qual_name[match.end() :]
                origin_full_name = f"{module_name}.{qual_name}"

            origin_name = QualifiedName.from_str(origin_full_name)

            for part in origin_name:
                if not part.isidentifier():
                    self.report_error(InvalidIdentifierError(part, path))
                    return None
            return origin_name

        def _parse_expression_str(self, expr_str: str) -> Value | InvalidExpression:
            strip_expr = expr_str.strip()
            try:
                ast.parse(strip_expr)
                print_safe = False
                try:
                    ast.literal_eval(strip_expr)
                    print_safe = True
                except (
                    ValueError,
                    TypeError,
                    SyntaxError,
                    MemoryError,
                    RecursionError,
                ):
                    pass
                return Value(strip_expr, is_print_safe=print_safe)
            except SyntaxError:
                self.report_error(InvalidExpressionError(strip_expr))
                return InvalidExpression(strip_expr)

        def parse_args_str(self, args_str: str) -> list[Argument]:
            split_args = self._split_args_str(args_str)
            if split_args is None:
                return [
                    Argument(name=Identifier("args"), variadic=True),
                    Argument(name=Identifier("kwargs"), kw_variadic=True),
                ]

            result: list[Argument] = []
            kw_only = False
            for arg_str, annotation_str, default_str in split_args:
                if arg_str.strip() == "/":
                    for arg in result:
                        arg.pos_only = True
                    continue
                if arg_str.strip() == "*":
                    kw_only = True
                    continue
                _arg_star_name_regex = re.compile(r"^\s*(?P<stars>\*{1,2})?" r"\s*(?P<name>\w+)\s*$")
                match = _arg_star_name_regex.match(arg_str)
                if match is None:
                    return [
                        Argument(name=Identifier("args"), variadic=True),
                        Argument(name=Identifier("kwargs"), kw_variadic=True),
                    ]
                name = Identifier(match.group("name"))

                variadic = False
                kw_variadic = False

                stars = match.group("stars")
                if stars == "*":
                    variadic = True
                elif stars == "**":
                    kw_variadic = True

                if annotation_str is not None:
                    annotation = self.parse_annotation_str(annotation_str)
                else:
                    annotation = None

                if default_str is not None:
                    default = self.parse_value_str(default_str)
                else:
                    default = None

                result.append(
                    Argument(
                        name=name,
                        default=default,
                        annotation=annotation,
                        variadic=variadic,
                        kw_variadic=kw_variadic,
                        kw_only=kw_only,
                    )
                )
            return result

        def parse_annotation_str(self, annotation_str: str) -> ResolvedType | InvalidExpression | Value:
            variants = self._split_type_union_str(annotation_str)
            if variants is None or len(variants) == 0:
                self.report_error(InvalidExpressionError(annotation_str))
                return InvalidExpression(annotation_str)
            if len(variants) == 1:
                return self.parse_type_str(variants[0])
            union_t = self.parse_annotation_str("typing.Union")
            assert isinstance(union_t, ResolvedType)
            return ResolvedType(
                name=union_t.name,
                parameters=[self.parse_type_str(variant) for variant in variants],
            )

        def parse_type_str(self, annotation_str: str) -> ResolvedType | InvalidExpression | Value:
            qname_regex = re.compile(r"^\s*(?P<qual_name>([_A-Za-z]\w*)?(\s*\.\s*[_A-Za-z]\w*)*)")
            annotation_str = annotation_str.strip()
            match = qname_regex.match(annotation_str)
            if match is None:
                return self.parse_value_str(annotation_str)
            qual_name = QualifiedName(Identifier(part) for part in match.group("qual_name").replace(" ", "").split("."))
            parameters_str = annotation_str[match.end("qual_name") :].strip()

            if len(parameters_str) == 0:
                parameters = None
            else:
                if parameters_str[0] != "[" or parameters_str[-1] != "]":
                    return self.parse_value_str(annotation_str)

                split_parameters = self._split_parameters_str(parameters_str[1:-1])
                if split_parameters is None:
                    return self.parse_value_str(annotation_str)

                parameters = [self.parse_annotation_str(param_str) for param_str in split_parameters]
            return ResolvedType(name=qual_name, parameters=parameters)

        def parse_function_docstring(self, func_name: Identifier, doc_lines: list[str]) -> list[Function]:
            if len(doc_lines) == 0:
                return []

            top_signature_regex = re.compile(rf"^{func_name}\((?P<args>.*)\)\s*(->\s*(?P<returns>.+))?$")

            match = top_signature_regex.match(doc_lines[0])
            if match is None:
                return []

            if len(doc_lines) < 2 or doc_lines[1] != "Overloaded function.":
                returns_str = match.group("returns")
                if returns_str is not None:
                    returns = self.parse_annotation_str(returns_str)
                else:
                    returns = None

                return [
                    Function(
                        name=func_name,
                        args=self.parse_args_str(match.group("args")),
                        doc=self._strip_empty_lines(doc_lines[1:]),
                        returns=returns,
                    )
                ]

            overload_signature_regex = re.compile(
                rf"^(\s*(?P<overload_number>\d+).\s*)"
                rf"{func_name}\((?P<args>.*)\)\s*->\s*(?P<returns>.+)$"
            )

            doc_start = 0
            _dummy = Function(Identifier(""))
            overloads = [_dummy]

            for i in range(2, len(doc_lines)):
                match = overload_signature_regex.match(doc_lines[i])
                if match:
                    if match.group("overload_number") != f"{len(overloads)}":
                        continue
                    overloads[-1].doc = self._strip_empty_lines(doc_lines[doc_start:i])
                    doc_start = i + 1
                    overloads.append(
                        Function(
                            name=func_name,
                            args=self.parse_args_str(match.group("args")),
                            returns=self.parse_annotation_str(match.group("returns")),
                            doc=None,
                            decorators=[Decorator(str(self.parse_annotation_str("typing.overload")))],
                        )
                    )

            overloads[-1].doc = self._strip_empty_lines(doc_lines[doc_start:])

            return overloads[1:]

        def _fixup_parsed_getters_or_setters(self, funcs: list[Function]) -> Function | None:
            if len(funcs) == 0:
                return None
            if len(funcs) > 1:
                raise RuntimeError("Multiple overloads in property's getters/setters are not supported")

            func = funcs[0]

            if (
                len(func.args) > 0
                and not func.args[0].variadic
                and not func.args[0].kw_variadic
                and func.args[0].default is None
            ):
                func.args[0].name = Identifier("self")
                func.args[0].annotation = None
            else:
                pass
            return func

        def _split_args_str(self, args_str: str) -> list[tuple[str, str | None, str | None]] | None:
            result = []

            closing = {"(": ")", "{": "}", "[": "]"}
            stack = []
            i = 0
            arg_begin = 0
            semicolon_pos: int | None = None
            eq_sign_pos: int | None = None

            def add_arg():
                nonlocal semicolon_pos
                nonlocal eq_sign_pos
                annotation = None
                default = None

                arg_end = i

                if eq_sign_pos is not None:
                    arg_end = eq_sign_pos
                    default = args_str[eq_sign_pos + 1 : i]

                if semicolon_pos is not None:
                    annotation = args_str[semicolon_pos + 1 : arg_end]
                    arg_end = semicolon_pos

                name = args_str[arg_begin:arg_end]
                result.append((name, annotation, default))
                semicolon_pos = None
                eq_sign_pos = None

            while i < len(args_str):
                c = args_str[i]
                if c in "\"'":
                    str_end = self._find_str_end(args_str, i)
                    if str_end is None:
                        return None
                    i = str_end
                elif c in closing:
                    stack.append(closing[c])
                elif len(stack) == 0:
                    if c == ",":
                        add_arg()
                        arg_begin = i + 1
                    elif c == ":" and semicolon_pos is None:
                        semicolon_pos = i
                    elif c == "=" and args_str[i : i + 2] != "==":
                        eq_sign_pos = i
                elif stack[-1] == c:
                    stack.pop()

                i += 1

            if len(stack) != 0:
                return None

            if len(args_str[arg_begin:i].strip()) != 0:
                add_arg()

            return result

        def _split_type_union_str(self, type_str: str) -> list[str] | None:
            return self._split_str(type_str, delim="|")

        def _split_parameters_str(self, param_str: str) -> list[str] | None:
            return self._split_str(param_str, delim=",")

        def _split_str(self, param_str: str, delim: str):
            result = []
            closing = {"(": ")", "{": "}", "[": "]"}
            stack = []
            i = 0
            arg_begin = 0

            def add_arg():
                arg_end = i
                param = param_str[arg_begin:arg_end]
                result.append(param)

            while i < len(param_str):
                c = param_str[i]
                if c in "\"'":
                    str_end = self._find_str_end(param_str, i)
                    if str_end is None:
                        return None
                    i = str_end
                elif c in closing:
                    stack.append(closing[c])
                elif len(stack) == 0:
                    if c == delim:
                        add_arg()
                        arg_begin = i + 1
                elif stack[-1] == c:
                    stack.pop()

                i += 1
            if len(stack) != 0:
                return None
            if param_str[arg_begin:i].strip() != 0:
                add_arg()
            return result

        def _find_str_end(self, s, start) -> int | None:
            for i in range(start + 1, len(s)):
                c = s[i]
                if c == "\\":
                    continue
                if c == s[start]:
                    return i
            return None

        def _strip_empty_lines(self, doc_lines: list[str]) -> Docstring | None:
            assert isinstance(doc_lines, list)
            start = 0
            for start in range(0, len(doc_lines)):
                if len(doc_lines[start].strip()) > 0:
                    break
            end = len(doc_lines) - 1
            for end in range(len(doc_lines) - 1, 0, -1):
                if len(doc_lines[end].strip()) > 0:
                    break
            result = "\n".join(doc_lines[start : end + 1])
            if len(result) == 0:
                return None
            return Docstring(result)

    class ExtractSignaturesFromPybind11Docstrings(IParser):
        _arg_star_name_regex = re.compile(r"^\s*(?P<stars>\*{1,2})?" r"\s*(?P<name>\w+)\s*$")

        def handle_function(self, path: QualifiedName, func: Any) -> list[Function]:
            result = super().handle_function(path, func)
            if len(result) == 1 and result[0].args == [
                Argument(name=Identifier("args"), variadic=True),
                Argument(name=Identifier("kwargs"), kw_variadic=True),
            ]:
                doc: str | None = func.__doc__
                func_name = Identifier(path[-1])

                if doc is not None:
                    parsed_funcs = self.parse_function_docstring(func_name, doc.splitlines())

                    if len(parsed_funcs) > 0:
                        return parsed_funcs
            return result

        def handle_property(self, path: QualifiedName, prop: Any) -> Property | None:
            result = Property(name=path[-1], modifier=None)

            fake_path = QualifiedName((*path, Identifier("")))

            if hasattr(prop, "fget") and prop.fget is not None:
                for func_path in [fake_path, path]:
                    result.getter = self._fixup_parsed_getters_or_setters(self.handle_function(func_path, prop.fget))
                    if result.getter is not None and result.getter.args != [
                        Argument(name=Identifier("args"), variadic=True),
                        Argument(name=Identifier("kwargs"), kw_variadic=True),
                    ]:
                        break

            if hasattr(prop, "fset") and prop.fset is not None:
                for func_path in [fake_path, path]:
                    result.setter = self._fixup_parsed_getters_or_setters(self.handle_function(func_path, prop.fset))
                    if result.setter is not None and result.setter.args != [
                        Argument(name=Identifier("args"), variadic=True),
                        Argument(name=Identifier("kwargs"), kw_variadic=True),
                    ]:
                        break
            if result.getter is None and result.setter is None:
                return None

            prop_doc = getattr(prop, "__doc__", None)
            if prop_doc is not None:
                result.doc = self._strip_empty_lines(prop_doc.splitlines())

            if (
                result.doc is not None
                and result.getter is not None
                and (
                    result.doc == result.getter.doc
                    or result.doc == self._strip_empty_lines((getattr(prop.fget, "__doc__", "") or "").splitlines())
                )
            ):
                result.doc = None

            return result

        def parse_args_str(self, args_str: str) -> list[Argument]:
            return BaseParser.parse_args_str(self, args_str)

        def parse_annotation_str(self, annotation_str: str) -> ResolvedType | InvalidExpression | Value:
            return BaseParser.parse_annotation_str(self, annotation_str)

        def parse_value_str(self, value: str) -> Value | InvalidExpression:
            return BaseParser.parse_value_str(self, value)

    class Parser(
        *error_handlers_top,
        FixMissing__future__AnnotationsImport,
        FixMissing__all__Attribute,
        FixMissingNoneHashFieldAnnotation,
        FixMissingImports,
        FilterTypingModuleAttributes,
        FixPEP585CollectionNames,
        FixTypingTypeNames,
        FixScipyTypeArguments,
        FixMissingFixedSizeImport,
        FixMissingEnumMembersAnnotation,
        OverridePrintSafeValues,
        *numpy_fixes,
        FixNumpyDtype,
        FixNumpyArrayFlags,
        FixCurrentModulePrefixInTypeNames,
        FixBuiltinTypes,
        RewritePybind11EnumValueRepr,
        FilterClassMembers,
        ReplaceReadWritePropertyWithField,
        FilterInvalidIdentifiers,
        FixValueReprRandomAddress,
        FixRedundantBuiltinsAnnotation,
        FilterPybindInternals,
        FilterPybind11ViewClasses,
        FixRedundantMethodsFromBuiltinObject,
        RemoveSelfAnnotation,
        FixPybind11EnumStrDoc,
        ExtractSignaturesFromPybind11Docstrings,
        ParserDispatchMixin,
        BaseParser,
        *error_handlers_bottom,
    ):
        pass

    parser = Parser()

    if args.enum_class_locations:
        parser.set_pybind11_enum_locations(dict(args.enum_class_locations))
    if args.ignore_invalid_identifiers is not None:
        parser.set_ignored_invalid_identifiers(args.ignore_invalid_identifiers)
    if args.ignore_invalid_expressions is not None:
        parser.set_ignored_invalid_expressions(args.ignore_invalid_expressions)
    if args.ignore_unresolved_names is not None:
        parser.set_ignored_unresolved_names(args.ignore_unresolved_names)
    if args.print_safe_value_reprs is not None:
        parser.set_print_safe_value_pattern(args.print_safe_value_reprs)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(name)s - [%(levelname)7s] %(message)s",
    )
    args = arg_parser().parse_args(argv, namespace=CLIArgs())

    parser = stub_parser_from_args(args)
    printer = Printer(invalid_expr_as_ellipses=not args.print_invalid_expressions_as_is)

    out_dir, sub_dir = to_output_and_subdir(
        output_dir=args.output_dir,
        module_name=args.module_name,
        root_suffix=args.root_suffix,
    )

    run(
        parser,
        printer,
        args.module_name,
        out_dir,
        sub_dir=sub_dir,
        dry_run=args.dry_run,
        writer=Writer(stub_ext=args.stub_extension),
    )


def to_output_and_subdir(output_dir: str, module_name: str, root_suffix: str | None) -> tuple[Path, Path | None]:
    out_dir = Path(output_dir)

    module_path = module_name.split(".")

    if root_suffix is None:
        return out_dir.joinpath(*module_path[:-1]), None
    else:
        module_path = [f"{module_path[0]}{root_suffix}", *module_path[1:]]
        if len(module_path) == 1:
            sub_dir = Path(module_path[-1])
        else:
            sub_dir = None
        return out_dir.joinpath(*module_path[:-1]), sub_dir


def run(
    parser: IParser,
    printer: Printer,
    module_name: str,
    out_dir: Path,
    sub_dir: Path | None,
    dry_run: bool,
    writer: Writer,
):
    module = parser.handle_module(QualifiedName.from_str(module_name), importlib.import_module(module_name))
    parser.finalize()

    if module is None:
        raise RuntimeError(f"Can't parse {module_name}")

    if dry_run:
        return

    out_dir.mkdir(exist_ok=True, parents=True)
    writer.write_module(module, printer, to=out_dir, sub_dir=sub_dir)


if __name__ == "__main__":
    main()
