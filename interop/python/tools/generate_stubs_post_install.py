#!/usr/bin/env python3
# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     class Email:
#         def __init__(self, v): self.v = v
#         def bind(self, f): return Email(f(self.v))
#         def get(self): return self.v
#     print(
#         Email("")
#         .bind(lambda s: s + "besian")
#         .bind(lambda s: s + "sejdiu")
#         .bind(lambda s: s + "@gmail.com")
#         .get()
#     )
#
"""
Stub generation script.

This script generates Python type stubs (.pyi files) for the lahuta C++ extension
modules directly in the local source tree. CMake will automatically copy them to
the install location alongside other Python files.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def _get_sanitizer_lib_path() -> str:
    if sys.platform != "darwin":
        return ""

    try:
        result = subprocess.run(
            ["clang", "--print-file-name=libclang_rt.tsan_osx_dynamic.dylib"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            lib_path = result.stdout.strip()
            if lib_path and "/" in lib_path:
                return lib_path
    except Exception:
        pass

    return ""


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: generate_stubs_post_install.py <install_prefix>", file=sys.stderr)
        return 1

    install_prefix = Path(sys.argv[1]).resolve()
    lahuta_pkg_dir = install_prefix / "lahuta"

    pkg_name_fmt = "[lahuta]"

    if not lahuta_pkg_dir.exists():
        print(f"{pkg_name_fmt} Package not found at {lahuta_pkg_dir}", file=sys.stderr)
        return 1

    # Determine local source directory (where stubs should be generated)
    pkg_root = Path(__file__).resolve().parents[1]
    local_lahuta_lib = pkg_root / "lahuta" / "lib"

    if not local_lahuta_lib.exists():
        print(f"{pkg_name_fmt} Local lib directory not found at {local_lahuta_lib}", file=sys.stderr)
        return 1

    print(f"{pkg_name_fmt} Generating stubs in local source tree", file=sys.stderr)

    stubgen_script = pkg_root / "external" / "pybind11_stubgen.py"

    if not stubgen_script.exists():
        print(f"{pkg_name_fmt} Stub generator not found at {stubgen_script}", file=sys.stderr)
        return 1

    env = os.environ.copy()
    prefix = str(install_prefix)
    env["PYTHONPATH"] = prefix + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    env.pop("PYTHONHOME", None)

    sanitizer_vars = ["DYLD_INSERT_LIBRARIES", "TSAN_OPTIONS", "LSAN_OPTIONS", "UBSAN_OPTIONS", "ASAN_OPTIONS"]
    for var in sanitizer_vars:
        if var in os.environ:
            env[var] = os.environ[var]

    if "DYLD_INSERT_LIBRARIES" not in env:
        lib_path = _get_sanitizer_lib_path()
        if lib_path:
            env["DYLD_INSERT_LIBRARIES"] = lib_path

    modules = ["lahuta.lib.lahuta"]

    mapping_module_path = lahuta_pkg_dir / "lib" / "mapping"
    if mapping_module_path.exists():
        modules.append("lahuta.lib.mapping")

    # pybind11-stubgen creates: <output-dir>/lahuta/lib/<module>/__init__.pyi
    # By setting output-dir to pkg_root (interop/python/), it creates the full path
    for module in modules:
        print(f"{pkg_name_fmt} Generating stubs for {module}", file=sys.stderr)

        try:
            result = subprocess.run(
                [
                    sys.executable,
                    str(stubgen_script),
                    module,
                    "--output-dir",
                    str(pkg_root),
                    "--root-suffix=",
                    "--ignore-all-errors",
                    "--numpy-array-use-type-var",
                ],
                env=env,
                cwd=str(install_prefix),
                capture_output=True,
                text=True,
                timeout=300,
            )

            if result.returncode != 0:
                print(f"{pkg_name_fmt} Stub generation failed for {module}:", file=sys.stderr)
                print(result.stderr, file=sys.stderr)
                return 1

        except Exception as e:
            print(f"{pkg_name_fmt} Error generating stubs for {module}: {e}", file=sys.stderr)
            return 1

    # Verify stubs were created
    module_name = "lahuta"
    generated_stub = local_lahuta_lib / module_name / "__init__.pyi"
    if generated_stub.exists():
        print(f"{pkg_name_fmt} Stubs generated at {local_lahuta_lib / module_name}", file=sys.stderr)
    else:
        print(f"{pkg_name_fmt} Warning: Expected stub not found at {generated_stub}", file=sys.stderr)

    print(f"{pkg_name_fmt} Stub generation completed", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
