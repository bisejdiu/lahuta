#!/usr/bin/env python3
"""Stub generation script."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


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

    print(f"{pkg_name_fmt} Generating stubs for package at {lahuta_pkg_dir}", file=sys.stderr)

    pkg_root = Path(__file__).resolve().parents[1]
    stubgen_script = pkg_root / "external" / "pybind11_stubgen.py"

    if not stubgen_script.exists():
        print(f"{pkg_name_fmt} Stub generator not found at {stubgen_script}", file=sys.stderr)
        return 1

    env = os.environ.copy()
    prefix = str(install_prefix)
    env["PYTHONPATH"] = prefix + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    env.pop("PYTHONHOME", None)

    modules = ["lahuta.lib.lahuta", "lahuta.lib.mapping"]
    for module in modules:
        print(f"{pkg_name_fmt} Generating stubs for {module}", file=sys.stderr)
        try:
            result = subprocess.run(
                [
                    sys.executable,
                    str(stubgen_script),
                    module,
                    "--output-dir",
                    str(install_prefix),
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

    print(f"{pkg_name_fmt} Stub generation completed successfully", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
