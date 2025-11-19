# ⚡ Lahuta ⚡
[![Github-CI][github-ci-core]][github-link]
[![Github-CI][github-ci-python]][github-link]
[![Conda][conda-badge]][conda-link]
[![Update Date][Anaconda-Update Badge]][update-link]

---

Structural Biology has entered a data-rich era.

**Lahuta** is a high-performance open-source software for scalable structural biology analysis. It creates chemically informative topological descriptors from diverse structural inputs and has been optimized for generative-AI model outputs (such as AlphaFold models). On a modern multi-core laptop, Lahuta processes the AlphaFold DB Swiss-Prot subset (~540,000 proteins) in minutes. In benchmarks, Lahuta is orders of magnitude faster than comparable tools.

Lahuta scales in both dataset and system size, handling hundreds of millions of structures as well as assemblies with tens of millions of atoms, and offers native support for molecular-dynamics simulations.

## Key Features

- **Scalability**: Processes hundreds of millions of structures and assemblies with tens of millions of atoms
- **Performance**: Sub- to low-millisecond computation times
- **Chemistry-aware topology**: Detailed topological representation with bond connectivity, bond orders, hybridization states, and protonation states
- **Format support**: PDB, PDBx/mmCIF, BinaryCIF, MMTF for structures; XTC and GRO for MD trajectories
- **MD simulation support**: Native, first-class support for analyzing molecular dynamics trajectories
- **Optimized for AI models**: Custom high-performance topology perception for AlphaFold models
- **Expressive selection language**: Query system with logical and arithmetic operators for precise substructure identification
- **Structural analysis**: Distance calculations, neighbor searching, native contact analysis
- **Extensive Contact Analysis**: Supports contact analysis based on most popular analysis tools (Arpeggio, MolStar, GetContacts)
- **Deep Python Integration**: Python is fully supported as a first-class interface and extendability layer
- **No HPC required**: Runs efficiently on standard laptops

## Design Principles

Lahuta development is guilded by the following three principles:
1. **High Performance**: All algorithms and data structures are (reasonably) optimized for speed.
2. **Scalability**: Designed to handle ultra large-scale datasets and systems with high efficiency.
3. **Foundational core library**: A composable, stable, high‑performance base that enables advanced structural analyses


## Motivation

Recent advances in structure prediction have resulted in an explosion of available structural data. The AlphaFold Database contains over 200 million models, and metagenomic predictions from ESMFold add over 600 million more—3-4 orders of magnitude larger than experimental archives. Meanwhile, new generators like BioEmu can synthesize vast conformational ensembles within hours. Existing tools in structural biology was not designed for this volume or heterogeneity, making chemistry-aware, ultra large-scale analyses either infeasible or dependent on prohibitively large compute resources.

Lahuta enables screening of millions of models and long MD ensembles on standard hardware, making it possible to discover novel conformational states, folds and fold families, systematic mapping of interfaces and ligandable pockets, and ensemble-level comparisons that were previously infeasible without dedicated compute resources.

## Installation

We recommend using `conda` for installation. Lahuta requires **Python 3.10** or higher.

```bash
conda install -c bisejdiu lahuta
```

Or `pip`:

```bash
pip install lahuta
```

## Build, Test, and Install (C++ Core, CLI, and Python bindings)

- Configure, build, and install from the repository root using Ninja:
  ```bash
  cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DLAHUTA_BUILD_PYTHON=ON -DLAHUTA_GENERATE_PY_STUBS=ON -DBUILD_TESTING=ON -DLAHUTA_BUILD_CLI=ON -DLAHUTA_BUILD_EXAMPLES=ON
  cmake --build build -j 8 && cmake --install build
  ```

- Run tests (via CTest):
  ```bash
  cd build && ctest --output-on-failure
  ```

- Notes and options:
  - `LAHUTA_BUILD_PYTHON=ON` builds the Python shared libraries and installs the Python package into the CMake install prefix. It does not install Python-level dependencies; use `pip` or `conda` to install those (see list in `interop/python/pyproject.toml`).
  - Switch between shared and static linkage for `lahuta_core` with:
    ```bash
    cmake -S . -B build -DLAHUTA_BUILD_SHARED_CORE=OFF
    ```
    The default (`ON`) produces a shared library for reuse by the Python bindings. Set to `OFF` for a static CLI.
  - Library-only builds (used by Python packaging) can disable the CLI:
    ```bash
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLAHUTA_BUILD_CLI=OFF
    ```
  - Optional test configurations:
    - `-DENABLE_ASAN=ON` - Enable AddressSanitizer and UndefinedBehaviorSanitizer
    - `-DENABLE_TSAN=ON` - Enable ThreadSanitizer
  - Run specific tests:
    ```bash
    cd build && ctest -R test_name_pattern --output-on-failure
    ```

## Project Structure Overview

The core of Lahuta consists of these modules (see under `core/src/`):

- `analysis` - Structural analysis routines including contact analysis, system-level properties, and topology perception.
- `bonds` - Bond perception algorithms, bond order assignment, and connectivity validation.
- `chemistry` - Chemical property calculations (formal charges, hydrophobicity, atom typing routines).
- `compute` - Computation abstraction with dependency management, pipeline execution, and result caching.
- `contacts` - Contact detection implementations from multiple methods (Arpeggio, MolStar, GetContacts).
- `db` - Database interface for fast storage and retrieval of structural data and zero-copy reads.
- `distances` - Distance calculations, neighbor search algorithms, and pairwise distance matrices.
- `entities` - Core entity representations (atoms, residues, contacts) and their associated views and iterators.
- `entities/search` - Entity query and retrieval with hit buffering.
- `md` - Molecular dynamics trajectory parsing (XTC, GRO formats) and frame-by-frame analysis support.
- `models` - Optimized system-level properties and topology perception for AI (currently AF2) models.
- `pipeline` - High-level pipeline framework for processing, parallel execution, backpressure system, and progress tracking.
- `selections` - Expressive selection language for querying atoms, residues, and substructures based on geometric, chemical, or topological criteria (WIP)
- `serialization` - Data serialization and deserialization.
- `sinks` - Data output handlers that write pipeline results to files, databases, or memory.
- `spatial` - Spatial indexing structures (cell lists, KD-trees) for scalable neighbor queries and contact searches.
- `topology` - Topology construction engine.

The project also includes:

- `cli/` - Command-line interface tools for structure analysis and database creation.
- `interop/python` - Python integration
- `core/tests/` - Comprehensive C++ test suite using Google Test framework.
- `core/examples/` - Example programs demonstrating C++ core library usage.

Python integration structure:
- `interop/python/src/` - Python bindings implementation (pybind11) with deep (often zero-copy) NumPy integration.
- `interop/python/lahuta/` - Python package providing high-level APIs, utilities, and type-safe interfaces to core functionality.
- `interop/python/examples/` - Example scripts demonstrating Python API usage. 
- `interop/python/tests/` - Python test suite (pytest) covering all public APIs and integration scenarios.
- `interop/python/benchmarks/` - Performance benchmarks comparing Lahuta against other tools.

## Documentation

See `interop/python/examples` and `interop/python/tests` for Python usage examples.

## Reporting Issues

Report issues in the [issues](https://github.com/bisejdiu/lahuta/issues) section.

## Contributing

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.


[github-ci-core]: https://github.com/bisejdiu/lahuta/actions/workflows/ci-core.yml/badge.svg?branch=main
[github-ci-python]: https://github.com/bisejdiu/lahuta/actions/workflows/ci-python.yml/badge.svg?branch=main
[github-link]: https://github.com/bisejdiu/lahuta
[conda-badge]: https://anaconda.org/bisejdiu/lahuta/badges/version.svg
[conda-link]: https://anaconda.org/bisejdiu/lahuta
[Anaconda-Update Badge]: https://anaconda.org/bisejdiu/lahuta/badges/latest_release_date.svg
[update-link]: https://anaconda.org/bisejdiu/lahuta
