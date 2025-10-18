# Lahuta Build Instructions (Developers)

This document describes how to build Lahuta from source for development purposes.

## Prerequisites

- CMake 3.15 or later
- C++17 compatible compiler
- Conda environment with required dependencies
- Python 3.8+ (for Python bindings)

## Environment Setup

```bash
# Activate the lahuta development environment
conda activate lahuta-cxx

# Set C++ include path for Boost dependencies
set -gx CPLUS_INCLUDE_PATH $CONDA_PREFIX/include
```

## Building the Core C++ Library and Executable

```bash
cd core/
cmake -S . -B build -G Ninja
cmake --build build -j 16
```

This produces:
- `build/lahuta` - Main executable
- `build/liblahuta_cpp.a` - Static C++ library
- `build/_lahuta.cpython-310-darwin.so` - Python extension (if pybind11 available)

## Python Bindings

The Python bindings are automatically built as part of the core build process if pybind11 is found.

### Automatic Installation
The CMake build automatically copies the Python extension to `interop/python/lahuta/lib/` for immediate use.

### Installing Python Bindings for Development
```bash
cd interop/python/
pip install -e .
```

### Testing Python Bindings
```bash
# From any directory:
python -c "import lahuta; print('Lahuta successfully imported')"

# Direct access to C++ classes:
python -c "from lahuta.lib._lahuta import FastNS_, NSResults_; print('Direct access works')"

# Or run the test suite:
cd interop/python/
python test_bindings.py
```

## Project Structure

```
lahuta/
├── core/                          # C++ core library and build system
│   ├── CMakeLists.txt             # Main build configuration
│   ├── src/                       # C++ source code
│   ├── external/                  # Third-party dependencies
│   └── build/                     # Build outputs
├── interop/                       # Language bindings
│   └── python/                    # Python-specific bindings
│       ├── src/                   # C++ binding source (pybind11)
│       ├── lahuta/                # Minimal Python package
│       │   ├── lib/               # Compiled extensions directory
│       │   │   └── _lahuta.*.so   # C++ extension module
│       │   └── __init__.py        # Import/verification functions
│       ├── pyproject.toml         # Python build configuration
│       └── test_bindings.py       # Binding tests
└── lahuta/                        # Main Python package (high-level API)
```

## Build Dependencies

The build automatically handles most dependencies via git submodules and CMake:
- Boost (graph library)
- pybind11 (Python bindings)
- Eigen3 (linear algebra)
- RDKit (chemistry toolkit)
- gemmi (crystallography)
- foldseek (structure search)
- distopia (distance calculations)

## Troubleshooting

### Missing Boost Headers
If you encounter missing `boost/iostreams/tee.hpp` or similar errors:
```bash
# Ensure C++ include path is set
set -gx CPLUS_INCLUDE_PATH $CONDA_PREFIX/include
```

### Python Import Errors
If Python can't find the `_lahuta` module:
1. Check that the build completed successfully
2. Verify the `.so` file exists in `interop/python/lahuta/lib/`
3. Run tests from the `interop/python/` directory

### Clean Build
```bash
cd core/
rm -rf build/
cmake -S . -B build -G Ninja
cmake --build build -j 16
```
