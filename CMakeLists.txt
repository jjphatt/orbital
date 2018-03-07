# Minimum CMake version.
cmake_minimum_required (VERSION 3.1.0)

# Look for the Polymec library.
set(MIN_POLYMEC_VERSION 2.0.0) # Minimum library version.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${POLYMEC_PREFIX}/share/polymec)
include(polymec)
if (NOT DEFINED POLYMEC_VERSION)
  message(FATAL_ERROR "Could not find polymec library. Please use the polymec=/path/to/polymec option.")
endif()
if (${POLYMEC_VERSION} STRLESS ${MIN_POLYMEC_VERSION})
  message(FATAL_ERROR "Orbital currently requires polymec v${MIN_POLYMEC_VERSION}, but you have installed v${POLYMEC_VERSION}.")
endif()
message("-- Found polymec library in ${POLYMEC_PREFIX} (v${POLYMEC_VERSION})")

set(CMAKE_C_COMPILER ${POLYMEC_C_COMPILER})
enable_language(C)

project (orbital)

# Build everything the way Polymec does.
set(CMAKE_C_FLAGS ${POLYMEC_C_FLAGS})
set(CMAKE_EXE_LINKER_FLAGS ${POLYMEC_EXE_LINKER_FLAGS})
set(BUILD_SHARED_LIBS ${POLYMEC_BUILD_SHARED_LIBS})
set(CMAKE_SHARED_LINKER_FLAGS ${POLYMEC_SHARED_LINKER_FLAGS})
set(CMAKE_BUILD_TYPE ${POLYMEC_BUILD_TYPE})
set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${POLYMEC_PREFIX}/lib")
enable_language(C)
set(ORBITAL_MAJOR_VERSION 0)
set(ORBITAL_MINOR_VERSION 1)
set(ORBITAL_PATCH_VERSION 0)
set(ORBITAL_VERSION "${ORBITAL_MAJOR_VERSION}.${ORBITAL_MINOR_VERSION}.${ORBITAL_PATCH_VERSION}")

# Install wherever Polymec lives.
set(CMAKE_INSTALL_PREFIX ${POLYMEC_PREFIX})

# Pass some additional flags.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DORBITAL_VERSION=\\\"${ORBITAL_VERSION}\\\"")

# Look for headers in these directories.
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_BINARY_DIR}/src")

enable_testing()

# The stuff.
add_subdirectory(src)

