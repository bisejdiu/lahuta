cmake_minimum_required(VERSION 3.15)

set(NATIVE_ARCH 1 CACHE BOOL "Assume native architecture for SIMD. Use one of the HAVE_* options or set CMAKE_CXX_FLAGS to the appropriate flags if you disable this.")

# project(lex2 CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -s")
# set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message("-- Compiler is clang(++)")
    set(CMAKE_COMPILER_IS_CLANG 1)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message("-- Compiler is GNU ")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.9.0")
        message(FATAL_ERROR "Insufficient gcc version")
    endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    message("-- Compiler is icc(++)")
    set(CMAKE_COMPILER_IS_ICC 1)
endif ()

# SIMD instruction sets support
set(MMSEQS_ARCH "")
if (HAVE_AVX2)
    if (CMAKE_COMPILER_IS_CLANG)
        set(MMSEQS_ARCH "${MMSEQS_ARCH} -mavx2 -mcx16")
    else ()
        set(MMSEQS_ARCH "${MMSEQS_ARCH} -mavx2 -mcx16 -Wa,-q")
    endif ()
    set(X64 1 CACHE INTERNAL "")
elseif (HAVE_SSE4_1)
    set(MMSEQS_ARCH "${MMSEQS_ARCH} -msse4.1 -mcx16")
    set(X64 1 CACHE INTERNAL "")
elseif (HAVE_SSE2)
    set(MMSEQS_ARCH "${MMSEQS_ARCH} -msse2")
    set(DISABLE_IPS4O 1)
    set(X64 1 CACHE INTERNAL "")
elseif (HAVE_POWER9)
    set(MMSEQS_ARCH "${MMSEQS_ARCH} -mcpu=power9 -mvsx")
    set(PPC64 1 CACHE INTERNAL "")
elseif (HAVE_POWER8)
    set(MMSEQS_ARCH "${MMSEQS_ARCH} -mcpu=power8 -mvsx")
    set(PPC64 1 CACHE INTERNAL "")
elseif (HAVE_ARM8)
    set(MMSEQS_ARCH "${MMSEQS_ARCH} -march=armv8-a+simd")
    set(ARM 1 CACHE INTERNAL "")
elseif (HAVE_S390X)
    set(MMSEQS_ARCH "${MMSEQS_ARCH} -march=z14 -mzarch -mzvector")
    set(ZARCH 1 CACHE INTERNAL "")
endif ()

if (NATIVE_ARCH AND (MMSEQS_ARCH STREQUAL ""))
    if (CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm.*|ARM.*|aarch64.*|AARCH64.*)")
        set(ARM 1 CACHE INTERNAL "")
    elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "PPC64*|ppc64*|powerpc64*")
        set(PPC64 1 CACHE INTERNAL "")
    elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "amd64|AMD64")
        set(X64 1 CACHE INTERNAL "")
    elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "x86|X86")
        set(X86 1 CACHE INTERNAL "")
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^s390")
        set(ZARCH 1 CACHE INTERNAL "")
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^sparc")
        set(SPARC 1 CACHE INTERNAL "")
    else ()
        message(WARNING "CPU without native SIMD instructions. Performance will be bad.")
    endif ()
    if (EMSCRIPTEN)
        set(MMSEQS_ARCH "-msimd128 -s WASM=1 -s ASSERTIONS=1")
    elseif (X86 OR X64)
        include(cmake/CheckSSEFeatures.cmake)
        if (NOT HAVE_SSE4_1_EXTENSIONS)
            if (NOT HAVE_SSE2_EXTENSIONS)
                message(FATAL_ERROR "At least SSE2 is needed to compile")
            endif ()
            message(WARNING "At least SSE4.1 is needed for best performance")
        endif ()
        # clang has a problem with march=native on travis
        if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.0.0")
            set(MMSEQS_ARCH "${SSE_FLAGS}")
        else()
            set(MMSEQS_ARCH "-march=native")
        endif()
    elseif (PPC64 OR ARM)
        set(MMSEQS_ARCH "-mcpu=native")
    elseif (ZARCH)
        set(MMSEQS_ARCH "-mcpu=native -mzvector")
    else ()
        set(MMSEQS_ARCH "-march=native")
    endif ()
endif ()

set(MMSEQS_ARCH ${MMSEQS_ARCH} CACHE INTERNAL "")

if (NOT (MMSEQS_ARCH STREQUAL ""))
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} ${MMSEQS_ARCH}")
endif ()

if (CYGWIN OR ARM OR PPC64)
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -D_GNU_SOURCE=1")
endif ()

if (CMAKE_COMPILER_IS_ICC)
    # default -fp-model results in inconsistent results in profile search
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -fp-model precise")
endif ()

# Apple specific features
if (APPLE)
    # macOS SDK started using _Atomic (in ucred.h) which g++ does not support
    # __APPLE_API_STRICT_CONFORMANCE makes sysctl.h not include apis like ucred.h
    # See: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=89864
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -D__APPLE_API_STRICT_CONFORMANCE")
endif ()

if (CMAKE_SYSTEM_NAME STREQUAL "FreeBSD")
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -D_WITH_GETLINE")
endif ()

if (CYGWIN)
    # default cygwin allocator (dlmalloc) locks on every allocation and destroys MT performance
    add_subdirectory(lib/nedmalloc)
    # IPS4O seems to deadlock on cygwin
    set(DISABLE_IPS4O 1)
endif()

if ((CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "5.0.0")
        OR (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.1.0"))
    # clang before v5 throws compile errors on ips4o
    set(DISABLE_IPS4O 1)
endif ()

if (PPC64 OR SPARC OR ZARCH)
    # FIXME: investigate why on ppc the regression seems to fail randomly
    set(DISABLE_IPS4O 1)
endif ()

set(MMSEQS_C_FLAGS "${MMSEQS_CXX_FLAGS}")
set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -std=c++1y")
# Compiler-specific features
if (CMAKE_COMPILER_IS_CLANG AND (NOT EMSCRIPTEN))
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++11")
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -stdlib=libc++")
endif ()

# include_directories(
#     ${CMAKE_CURRENT_LIST_DIR}
#     ${CMAKE_CURRENT_LIST_DIR}/3di
#     ${CMAKE_CURRENT_LIST_DIR}/alp
#     ${CMAKE_CURRENT_LIST_DIR}/cacode
#     ${CMAKE_CURRENT_LIST_DIR}/gemmi
#     ${CMAKE_CURRENT_LIST_DIR}/generated
#     ${CMAKE_CURRENT_LIST_DIR}/simd
#     ${CMAKE_CURRENT_LIST_DIR}/tmalign
# )

set(ENV{CARGO_NET_OFFLINE} true)
add_subdirectory(corrosion)
# don't try to link to gcc_s, its not needed anyway
list(REMOVE_ITEM Rust_CARGO_TARGET_LINK_NATIVE_LIBS "gcc_s")

if(Rust_VERSION VERSION_LESS "1.68.0" AND NOT IGNORE_RUST_VERSION)
    message(FATAL_ERROR "Rust >= 1.68.0 is required")
endif()

set(RUST_FEATURE "")
if(HAVE_AVX2)
    set(RUST_FEATURE simd_avx2)
elseif(HAVE_SSE4_1 OR HAVE_SSE2)
    set(RUST_FEATURE simd_sse2)
elseif(HAVE_ARM8)
    set(RUST_FEATURE simd_neon)
endif()

if (NATIVE_ARCH AND (RUST_FEATURE STREQUAL ""))
    if(EMSCRIPTEN)
        set(RUST_FEATURE simd_wasm)
    elseif(ARM)
        set(RUST_FEATURE simd_neon)
    elseif(X86 OR X64)
        if (HAVE_AVX2_EXTENSIONS)
            set(RUST_FEATURE simd_avx2)
        elseif (HAVE_SSE2_EXTENSIONS)
            set(RUST_FEATURE simd_sse2)
        endif()
    endif()
endif()

if(RUST_FEATURE STREQUAL "")
    message(FATAL_ERROR "No SIMD support for this architecture")
endif()
message("-- Rust Feature: ${RUST_FEATURE}")
corrosion_import_crate(
    MANIFEST_PATH block-aligner/c/Cargo.toml
    FEATURES "${RUST_FEATURE}"
    CRATE_TYPES staticlib
)
include_directories(
    ${CMAKE_CURRENT_LIST_DIR}/block-aligner/c
)

add_library(foldseek SHARED
    BaseMatrix.cpp
    CacheFriendlyOperations.cpp
    Debug.cpp
    EvalueNeuralNet.cpp
    Indexer.cpp
    KmerGenerator.cpp
    LDDT.cpp
    QueryMatcher.cpp
    Sequence.cpp
    SequenceLookup.cpp
    StructureSmithWaterman.cpp
    SubstitutionMatrix.cpp
    tantan.cpp
    TMaligner.cpp
    UngappedAlignment.cpp
    Util.cpp

    ${CMAKE_CURRENT_LIST_DIR}/cacode/lambda_calculator.cpp
    ${CMAKE_CURRENT_LIST_DIR}/cacode/lubksb.cpp
    ${CMAKE_CURRENT_LIST_DIR}/cacode/ludcmp.cpp
    ${CMAKE_CURRENT_LIST_DIR}/cacode/nrutil.cpp
    ${CMAKE_CURRENT_LIST_DIR}/kerasify/keras_model.cpp
    ${CMAKE_CURRENT_LIST_DIR}/3di/structureto3di.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tmalign/TMalign.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tmalign/affineneedlemanwunsch.cpp
)

target_include_directories(foldseek PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
    ${CMAKE_CURRENT_LIST_DIR}/3di
    ${CMAKE_CURRENT_LIST_DIR}/alp
    ${CMAKE_CURRENT_LIST_DIR}/cacode
    ${CMAKE_CURRENT_LIST_DIR}/block-aligner/c
    ${CMAKE_CURRENT_LIST_DIR}/generated
    ${CMAKE_CURRENT_LIST_DIR}/simd
    ${CMAKE_CURRENT_LIST_DIR}/tmalign
    $ENV{CONDA_PREFIX}/include
)

set_source_files_properties(
    ${CMAKE_CURRENT_LIST_DIR}/tmalign/TMalign.cpp
    PROPERTIES COMPILE_FLAGS -Wno-absolute-value
)

target_link_libraries(foldseek PRIVATE gemmi_cpp)

# add_executable(lex2 ${SOURCES})
# target_link_libraries(lex2 PRIVATE block_aligner_c)

# find_package(OpenMP REQUIRED)
# if(OpenMP_CXX_FOUND)
#     target_link_libraries(lex2 PRIVATE OpenMP::OpenMP_CXX)
# endif()
target_link_libraries(foldseek PUBLIC lahuta_shared_headers)
target_link_libraries(foldseek PRIVATE block_aligner_c)
set_property(TARGET foldseek PROPERTY POSITION_INDEPENDENT_CODE ON)

