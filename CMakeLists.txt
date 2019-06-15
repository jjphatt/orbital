# Minimum CMake version.
cmake_minimum_required (VERSION 3.10.0)

# Look for the scasm library.
set(MIN_SCASM_VERSION 3.0.0) # Minimum library version.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${SCASM_PREFIX}/share/scasm)
include(scasm)
if (NOT DEFINED SCASM_VERSION)
  message(FATAL_ERROR "Could not find the scasm library. Please set SCASM in config.sh.")
endif()
if (${SCASM_VERSION} STRLESS ${MIN_SCASM_VERSION})
  message(FATAL_ERROR "Orbital currently requires scasm v${MIN_SCASM_VERSION}, but you have installed v${SCASM_VERSION}.")
endif()
message("-- Found scasm library in ${SCASM_PREFIX} (v${SCASM_VERSION})")

set(CMAKE_C_COMPILER ${SCASM_C_COMPILER})
enable_language(C)

project (orbital)

# Build everything the way SCASM does.
set(CMAKE_C_FLAGS ${SCASM_C_FLAGS})
set(CMAKE_EXE_LINKER_FLAGS ${SCASM_EXE_LINKER_FLAGS})
set(BUILD_SHARED_LIBS ${SCASM_BUILD_SHARED_LIBS})
set(CMAKE_SHARED_LINKER_FLAGS ${SCASM_SHARED_LINKER_FLAGS})
set(CMAKE_BUILD_TYPE ${SCASM_BUILD_TYPE})
set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${SCASM_PREFIX}/lib")
enable_language(C)
set(ORBITAL_MAJOR_VERSION 0)
set(ORBITAL_MINOR_VERSION 1)
set(ORBITAL_PATCH_VERSION 0)
set(ORBITAL_VERSION "${ORBITAL_MAJOR_VERSION}.${ORBITAL_MINOR_VERSION}.${ORBITAL_PATCH_VERSION}")

# Install wherever scasm lives.
set(CMAKE_INSTALL_PREFIX ${SCASM_PREFIX})

# Pass some additional flags.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DORBITAL_VERSION=\\\"${ORBITAL_VERSION}\\\"")

# Look for headers in these directories.
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_BINARY_DIR}/src")

enable_testing()

# The stuff.
add_subdirectory(src)

