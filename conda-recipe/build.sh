#!/bin/bash
set -ex

BUILD_DIR="build-conda"

# Build with CMake directly to include both Python bindings and CLI
cmake ${CMAKE_ARGS} -S . -B "${BUILD_DIR}" -G Ninja \
  -DBUILD_TESTING=OFF \
  -DLAHUTA_BUILD_CLI=ON \
  -DLAHUTA_BUILD_SHARED_CORE=ON \
  -DLAHUTA_BUILD_PYTHON=ON \
  -DLAHUTA_GENERATE_PY_STUBS=OFF \
  -DLAHUTA_VENDOR_BOOST=ON \
  -DLAHUTA_USE_VENDORED_ZLIB_NG=ON \
  -DPython3_EXECUTABLE="${PYTHON}"

cmake --build "${BUILD_DIR}" --parallel "${CPU_COUNT}"

# For conda, install the Python package directly into site-packages from the
# existing CMake build tree. This avoids a second full scikit-build rebuild.
cmake --install "${BUILD_DIR}" --component python --prefix "${SP_DIR}"
cmake --install "${BUILD_DIR}" --component cli
