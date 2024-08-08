cmake_minimum_required(VERSION 3.15...3.22)
project(lahuta VERSION 1.0.0 LANGUAGES CXX C)

include(GNUInstallDirs)  # for CMAKE_INSTALL_LIBDIR, etc

if (DEFINED SKBUILD)  # building with scikit-build-core (pip install)
  # we set wheel.install-dir="/data" in pyproject.toml,
  # so here we need to only set paths to /platlib and /scripts
  set(PYTHON_INSTALL_DIR "${SKBUILD_PLATLIB_DIR}")
endif()

option(BUILD_SHARED_LIBS "Build shared library" ON)
option(USE_PYTHON "Build Python bindings" ON)
# option(USE_ZLIB_NG "Use zlib-ng instead of zlib" OFF)
# option(FETCH_ZLIB_NG "Download, build (statically) and use zlib-ng" ON)
# option(INTERNAL_ZLIB "Use subset of zlib distributed with gemmi" OFF)
option(GENERATE_STUBS "Generate Python type stubs" OFF)
option(STANDALONE_PYTHON_MODULE "Avoid linking Python module to liblahuta_cpp DLL" ON)

if (DEFINED Python_EXECUTABLE)
  set(USE_PYTHON ON CACHE BOOL "Build Python bindings" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
include(CheckIncludeFile)

if (DEFINED ENV{CXXFLAGS})
  set(USING_ENV_CXXFLAGS ON CACHE BOOL "" FORCE)
endif()

# If CMAKE_BUILD_TYPE is not set, and CXXFLAGS is not set in the environment,
# set it to Release. This is to avoid the default Debug build type.
if (NOT CMAKE_BUILD_TYPE AND NOT USING_ENV_CXXFLAGS)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# build using -g3 -Og if Debug, -O3 otherwise
# if (CMAKE_BUILD_TYPE STREQUAL "Debug")
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -Og")
# else()
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
# endif()


# if(NOT DEFINED RDKit_ROOT)
#     if(DEFINED ENV{CONDA_PREFIX})
#         set(RDKit_ROOT $ENV{CONDA_PREFIX})
#     else()
#         message(FATAL_ERROR "RDKit_ROOT is not defined. Please set it to the RDKit installation directory.")
#     endif()
# endif()

# list(APPEND CMAKE_PREFIX_PATH ${RDKit_ROOT}) 

# find_package(RDKit REQUIRED)

add_subdirectory(external/gemmi)
add_subdirectory(external/rdkit)

add_library(lahuta_cpp 
            src/bonds.cpp src/bond_order.cpp src/convert.cpp src/nsgrid.cpp 
            src/ob/kekulize.cpp src/ob/clean_mol.cpp external/ob/bitvec.cpp

            # external/rdkit/GraphMol/atomic_data.cpp
            # external/rdkit/GraphMol/PeriodicTable.cpp
            #
            # external/rdkit/DataStructs/base64.cpp external/rdkit/DataStructs/BitOps.cpp
            # external/rdkit/DataStructs/BitVect.cpp external/rdkit/DataStructs/ExplicitBitVect.cpp
            # external/rdkit/DataStructs/SparseBitVect.cpp 
            # external/rdkit/Geometry/point.cpp
            # external/rdkit/Geometry/Transform3D.cpp
            # external/rdkit/GraphMol/Atom.cpp 
            #
            # external/rdkit/GraphMol/AtomIterators.cpp external/rdkit/GraphMol/Atropisomers.cpp
            # external/rdkit/GraphMol/Bond.cpp external/rdkit/GraphMol/BondIterators.cpp
            # external/rdkit/GraphMol/Canon.cpp external/rdkit/GraphMol/Chirality.cpp
            # external/rdkit/GraphMol/Conformer.cpp
            # external/rdkit/GraphMol/Kekulize.cpp
            # external/rdkit/GraphMol/NontetrahedralStereo.cpp
            # external/rdkit/GraphMol/AddHs.cpp
            # external/rdkit/GraphMol/ConjugHybrid.cpp
            # external/rdkit/GraphMol/Aromaticity.cpp
            # external/rdkit/GraphMol/FindStereo.cpp
            # external/rdkit/GraphMol/Matrices.cpp
            # external/rdkit/GraphMol/AdjustQuery.cpp
            # external/rdkit/GraphMol/WedgeBonds.cpp
            # external/rdkit/GraphMol/FindRings.cpp
            #
            # external/rdkit/GraphMol/GraphMol.cpp external/rdkit/GraphMol/MolOps.cpp
            # external/rdkit/GraphMol/MonomerInfo.cpp external/rdkit/GraphMol/new_canon.cpp
            # external/rdkit/GraphMol/QueryAtom.cpp
            # external/rdkit/GraphMol/QueryBond.cpp external/rdkit/GraphMol/QueryOps.cpp
            # external/rdkit/GraphMol/Resonance.cpp external/rdkit/GraphMol/RingInfo.cpp
            # external/rdkit/GraphMol/ROMol.cpp external/rdkit/GraphMol/RWMol.cpp
            # external/rdkit/GraphMol/StereoGroup.cpp external/rdkit/GraphMol/SubstanceGroup.cpp
            #
            # external/rdkit/GraphMol/CIPLabeler/CIPLabeler.cpp external/rdkit/GraphMol/CIPLabeler/CIPMol.cpp
            # external/rdkit/GraphMol/CIPLabeler/Digraph.cpp external/rdkit/GraphMol/CIPLabeler/Edge.cpp
            # external/rdkit/GraphMol/CIPLabeler/Mancude.cpp external/rdkit/GraphMol/CIPLabeler/Node.cpp
            # external/rdkit/GraphMol/CIPLabeler/Sort.cpp 
            # external/rdkit/GraphMol/CIPLabeler/configs/AtropisomerBond.cpp
            # external/rdkit/GraphMol/CIPLabeler/configs/Configuration.cpp
            # external/rdkit/GraphMol/CIPLabeler/configs/Sp2Bond.cpp
            # external/rdkit/GraphMol/CIPLabeler/configs/Tetrahedral.cpp
            #
            # external/rdkit/GraphMol/CIPLabeler/rules/Rule1a.cpp
            # external/rdkit/GraphMol/CIPLabeler/rules/Rule1b.cpp
            # external/rdkit/GraphMol/CIPLabeler/rules/Rule2.cpp
            # external/rdkit/GraphMol/CIPLabeler/rules/Rule3.cpp
            # external/rdkit/GraphMol/CIPLabeler/rules/Rule4a.cpp
            # external/rdkit/GraphMol/CIPLabeler/rules/Rule4b.cpp
            # external/rdkit/GraphMol/CIPLabeler/rules/Rule4c.cpp
            # external/rdkit/GraphMol/CIPLabeler/rules/Rule5.cpp
            # external/rdkit/GraphMol/CIPLabeler/rules/Rule5New.cpp
            # external/rdkit/GraphMol/CIPLabeler/rules/Rule6.cpp
            # external/rdkit/GraphMol/CIPLabeler/rules/SequenceRule.cpp
            #
            # external/rdkit/GraphMol/GenericGroups/GenericGroups.cpp
            # external/rdkit/GraphMol/GraphMol.cpp external/rdkit/GraphMol/CIPLabeler/Mancude.cpp
            # external/rdkit/GraphMol/FileParsers/MolFileStereochem.cpp
            # external/rdkit/GraphMol/MolOps.cpp external/rdkit/GraphMol/MonomerInfo.cpp
            # external/rdkit/GraphMol/new_canon.cpp external/rdkit/GraphMol/CIPLabeler/Node.cpp
            # external/rdkit/GraphMol/QueryAtom.cpp
            # external/rdkit/GraphMol/QueryBond.cpp external/rdkit/GraphMol/QueryOps.cpp
            # external/rdkit/GraphMol/Resonance.cpp external/rdkit/GraphMol/RingInfo.cpp
            # external/rdkit/GraphMol/ROMol.cpp external/rdkit/GraphMol/RWMol.cpp
            #
            # external/rdkit/GraphMol/SmilesParse/SmilesParse.cpp
            # external/rdkit/GraphMol/SmilesParse/SmilesParseOps.cpp
            # external/rdkit/GraphMol/SmilesParse/CXSmilesOps.cpp
            # external/rdkit/GraphMol/SmilesParse/SmilesWrite.cpp
            #
            # external/rdkit/GraphMol/SmilesParse/smarts.tab.cpp
            # external/rdkit/GraphMol/SmilesParse/smiles.tab.cpp
            #
            # external/rdkit/GraphMol/SmilesParse/lex.yysmiles.cpp
            # external/rdkit/GraphMol/SmilesParse/lex.yysmarts.cpp
            #
            # external/rdkit/GraphMol/CIPLabeler/Sort.cpp
            # external/rdkit/GraphMol/StereoGroup.cpp external/rdkit/GraphMol/SubstanceGroup.cpp
            # external/rdkit/GraphMol/Substruct/SubstructMatch.cpp
            # external/rdkit/GraphMol/Substruct/SubstructUtils.cpp
            # external/rdkit/RDGeneral/Invariant.cpp external/rdkit/RDGeneral/LocaleSwitcher.cpp
            # external/rdkit/RDGeneral/RDLog.cpp external/rdkit/RDGeneral/types.cpp
            # external/rdkit/RDGeneral/utils.cpp

)

set_property(TARGET lahuta_cpp PROPERTY POSITION_INDEPENDENT_CODE ON)
target_compile_definitions(lahuta_cpp PUBLIC LAHUTA_BUILD)
 
target_include_directories(lahuta_cpp PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/external/ob
  # ${CMAKE_CURRENT_SOURCE_DIR}/external/rdkit
  # ${CMAKE_CURRENT_SOURCE_DIR}/external/rdkit/DataStructs
  # ${CMAKE_CURRENT_SOURCE_DIR}/external/rdkit/Geometry
  # ${CMAKE_CURRENT_SOURCE_DIR}/external/rdkit/GraphMol
  # ${CMAKE_CURRENT_SOURCE_DIR}/external/rdkit/GraphMol/CIPLabeler
  # ${CMAKE_CURRENT_SOURCE_DIR}/external/rdkit/GraphMol/GenericGroups
  # ${CMAKE_CURRENT_SOURCE_DIR}/external/rdkit/GraphMol/FileParsers
  # ${CMAKE_CURRENT_SOURCE_DIR}/external/rdkit/GraphMol/SmilesParse
  # ${CMAKE_CURRENT_SOURCE_DIR}/external/rdkit/GraphMol/Substruct
  # ${CMAKE_CURRENT_SOURCE_DIR}/external/rdkit/GraphMol/CIPLabeler/rules
  # ${CMAKE_CURRENT_SOURCE_DIR}/external/rdkit/GraphMol/CIPLabeler/configs
  # ${CMAKE_CURRENT_SOURCE_DIR}/external/rdkit/RDGeneral
)

# target_link_libraries(lahuta_cpp PUBLIC
#     $ENV{CONDA_PREFIX}/lib/libRDKitGraphMol.so
#     $ENV{CONDA_PREFIX}/lib/libRDKitRDGeneral.so
#     $ENV{CONDA_PREFIX}/lib/libRDKitRDGeometryLib.so
#     $ENV{CONDA_PREFIX}/lib/libRDKitSubstructMatch.so
#     $ENV{CONDA_PREFIX}/lib/libRDKitSmilesParse.so
# )

add_executable(lex main.cpp)
target_link_libraries(lex PUBLIC lahuta_cpp gemmi_cpp rdkit)

# install(TARGETS lahuta_cpp)
# install(TARGETS lahuta_cpp EXPORT lahuta_cpp
#         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
#         INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
#
# install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
#
# install(EXPORT lahuta_cpp
#         FILE lahuta_cpp-config.cmake
#         NAMESPACE lahuta_cpp::
#         DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lahuta_cpp)
#
# export(TARGETS lahuta_cpp NAMESPACE lahuta_cpp:: FILE lahuta_cpp-config.cmake)



### Python bindings ###
if (USE_PYTHON)
  message(STATUS "Building Python module.")
  # CMake >=3.18 has subcomponent Development.Module, scikit-build-core also has it
  if (${CMAKE_VERSION} VERSION_LESS 3.18 AND NOT SKBUILD)
    find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development)
  else()
    find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development.Module)
  endif()

  # use pybind11-config (if available) to determine pybind11_DIR
  execute_process(COMMAND pybind11-config --cmakedir OUTPUT_VARIABLE pybind11_DIR)
  string(STRIP "${pybind11_DIR}" pybind11_DIR)
  find_package(pybind11 2.6 CONFIG REQUIRED)
  message(STATUS "Found pybind11 ${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")

  pybind11_add_module(lahuta_py python/lahuta.cpp)
  # target_link_libraries(lahuta_py PRIVATE lahuta_cpp gemmi_cpp)

  if (STANDALONE_PYTHON_MODULE)
    target_sources(lahuta_py PRIVATE 
    $<TARGET_OBJECTS:lahuta_cpp>
    $<TARGET_OBJECTS:gemmi_cpp>
    )
    get_target_property(_lahuta_cpp_libs lahuta_cpp LINK_LIBRARIES)
    get_target_property(_gemmi_cpp_libs gemmi_cpp LINK_LIBRARIES)
    set_property(TARGET lahuta_py APPEND PROPERTY LINK_LIBRARIES ${_lahuta_cpp_libs})
    set_property(TARGET lahuta_py APPEND PROPERTY LINK_LIBRARIES ${_gemmi_cpp_libs})
  else()
    target_link_libraries(lahuta_py PRIVATE lahuta_cpp gemmi_cpp)
  endif()

  set_property(TARGET lahuta_py PROPERTY OUTPUT_NAME _lahuta)
  if(CMAKE_CXX_FLAGS MATCHES "-Wshadow")
    target_compile_options(lahuta_py PRIVATE "-Wno-shadow")
  endif()

  if (GENERATE_STUBS AND CMAKE_CROSSCOMPILING)
    message(WARNING "Stubs cannot be generated when cross-compiling - skipping.")
  elseif (GENERATE_STUBS)
    set(sep ":")
    if (WIN32)
      set(sep ";")
    endif()
    add_custom_command(TARGET lahuta_py POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E env
         "PYTHONPATH=$<TARGET_FILE_DIR:lahuta_py>${sep}$ENV{PYTHONPATH}"
       ${Python_EXECUTABLE} -m pybind11_stubgen _lahuta
     COMMAND ${CMAKE_COMMAND} -E touch stubs/py.typed
     COMMENT "Generate type stubs with pybind11-stubgen"
     BYPRODUCTS stubs/__init__.pyi stubs/_lahuta.pyi stubs/py.typed
     VERBATIM)
  endif()
else()
  message(STATUS "Skipping Python module. Add -D USE_PYTHON=1 to build it.")
endif()


# if (USE_PYTHON)
#   if (DEFINED PYTHON_INSTALL_DIR)
#     message(STATUS "Install directory for Python module: ${PYTHON_INSTALL_DIR}")
#     set(Python_SITEARCH "${PYTHON_INSTALL_DIR}")
#   endif()
#   file(TO_CMAKE_PATH "${Python_SITEARCH}" Python_SITEARCH)
#   # Using Python_SITEARCH/SITELIB is not good, because they are absolute
#   # and don't respect CMAKE_INSTALL_PREFIX.
#   # https://discourse.cmake.org/t/findpython3-how-to-specify-local-installation-directory-for-python-module/3580/5
#   if (BUILD_SHARED_LIBS AND NOT STANDALONE_PYTHON_MODULE)
#     if (DEFINED SKBUILD)
#       # When building via pip, Python_SITEARCH is /tmp.../wheel/platlib,
#       # different from the actual install path. Best guess:
#       set(rpath_rel "../../../${CMAKE_INSTALL_LIBDIR}")
#     else()
#       file(RELATIVE_PATH rpath_rel
#            "${Python_SITEARCH}" "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
#     endif()
#     set_target_properties(lahuta_py PROPERTIES INSTALL_RPATH "${origin}/${rpath_rel}")
#   endif()
#   install(TARGETS lahuta_py DESTINATION "${Python_SITEARCH}" COMPONENT py)
#   # install(DIRECTORY examples DESTINATION "${Python_SITEARCH}/_lahuta" COMPONENT py
#   #         FILES_MATCHING PATTERN "*.py"
#   #         PATTERN "[._]*" EXCLUDE)
#   if (GENERATE_STUBS AND NOT CMAKE_CROSSCOMPILING)
#     install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/stubs/_lahuta"
#             DESTINATION "${Python_SITEARCH}" COMPONENT py)
#   endif()
# endif()
#
