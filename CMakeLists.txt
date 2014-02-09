# Author: petter.strandmark@gmail.com (Petter Strandmark)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0)
PROJECT(curve_extraction C CXX)

OPTION(BUILD_TESTING "Enable tests" ON)
OPTION(BUILD_EXAMPLES "Build examples" ON)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Default locations to search for on various platforms.
LIST(APPEND SEARCH_LIBS /usr/lib)
LIST(APPEND SEARCH_LIBS /usr/local/lib)
LIST(APPEND SEARCH_LIBS /usr/local/homebrew/lib) # Mac OS X
LIST(APPEND SEARCH_LIBS /opt/local/lib)
LIST(APPEND SEARCH_HEADERS "C:/Program files/spii/lib")

LIST(APPEND SEARCH_HEADERS /usr/include)
LIST(APPEND SEARCH_HEADERS /usr/local/include)
LIST(APPEND SEARCH_HEADERS /usr/local/homebrew/include) # Mac OS X
LIST(APPEND SEARCH_HEADERS /opt/local/include)
LIST(APPEND SEARCH_HEADERS "C:/Program Files/SPII/include")

ENABLE_TESTING()

# Multithreading using OpenMP
OPTION(OPENMP
       "Enable multi-threading (requires OpenMP)"
       ON)

IF (${OPENMP})
  FIND_PACKAGE(OpenMP)
  IF(${OPENMP_FOUND})
    MESSAGE("-- Found OpenMP.")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    ADD_DEFINITIONS(-DUSE_OPENMP)
  ELSE(${OPENMP_FOUND})
    MESSAGE("-- Can't find OpenMP. Continuing without it.")
  ENDIF(${OPENMP_FOUND})
ENDIF (${OPENMP})

FIND_PATH(SPII_INCLUDE NAMES spii-thirdparty/badiff.h PATHS ${SEARCH_HEADERS})
IF (NOT EXISTS ${SPII_INCLUDE})
  MESSAGE(FATAL_ERROR
          "Can't find Spii include directory. You should specify the location manually or install it.")
ELSE (NOT EXISTS ${SPII_INCLUDE})
  MESSAGE("-- Found Spii include directory.")
ENDIF (NOT EXISTS ${SPII_INCLUDE})

INCLUDE_DIRECTORIES(
  include
  ${SPII_INCLUDE}
  thirdparty/Catch
  thirdparty/Eigen
  thirdparty/gheap
  thirdparty/pgm
  )

FILE(GLOB CURVE_EXTRACTION_HEADERS ${CMAKE_SOURCE_DIR}/include/curve_extraction/*.h)
INSTALL(FILES ${CURVE_EXTRACTION_HEADERS} DESTINATION include/curve_extraction)

# Library dependencies
SET (CURVE_EXTRACTION_LIBRARY_DEPENDENCIES)

# OpenMP needs to be included as a library on some platforms.
IF (${OPENMP_FOUND})
  IF (NOT MSVC)
    LIST(APPEND CURVE_EXTRACTION_LIBRARY_DEPENDENCIES gomp)
  ENDIF (NOT MSVC)
ENDIF (${OPENMP_FOUND})


# Change the default build type from Debug to Release, while still
# supporting overriding the build type.
#
# The CACHE STRING logic here and elsewhere is needed to force CMake
# to pay attention to the value of these variables.
IF (NOT CMAKE_BUILD_TYPE)
  MESSAGE("-- No build type specified; defaulting to CMAKE_BUILD_TYPE=Release.")
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

# C++11 support.
include(EnableCPP11.cmake)

IF (CMAKE_COMPILER_IS_GNUCXX)
  # GCC is not strict enough by default, so enable most of the warnings.
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic -Wall -Wextra -Wno-unknown-pragmas -Wno-sign-compare -Wno-unused -Wno-unused-parameter")

  IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
    MESSAGE("-- Debug mode enabled for Gcc; adding support for Gcov.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
  ENDIF (CMAKE_BUILD_TYPE STREQUAL "Debug")
ENDIF (CMAKE_COMPILER_IS_GNUCXX)

#
# Clang settings
#
IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # Enable warnings for Clang.
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic -Wall -Wextra -Wno-sign-compare -Wno-unused-parameter")
ENDIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

# We're using std::functions with 6 arguments.
IF (MSVC)
  # Disable deprecation warning for standard functions.
  ADD_DEFINITIONS("/wd4996")
  # Loss of data
  ADD_DEFINITIONS("/wd4244")
  ADD_DEFINITIONS("-D_VARIADIC_MAX=6")
  ADD_DEFINITIONS("-DNOMINMAX=1")
ENDIF (MSVC)

IF (CMAKE_BUILD_TYPE STREQUAL "Release")
  IF (CMAKE_COMPILER_IS_GNUCXX)
      MESSAGE("-- Optimizing for release mode.")
      # Linux
      IF (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        MESSAGE("-- Linux")
        ADD_DEFINITIONS("-march=native -mtune=native")
      ENDIF (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
      # Cygwin
      IF (${CMAKE_SYSTEM_NAME} STREQUAL "CYGWIN")
        MESSAGE("-- Cygwin")
        ADD_DEFINITIONS("-march=native -mtune=native")
      ENDIF (${CMAKE_SYSTEM_NAME} STREQUAL "CYGWIN")
      # Mac OS X
      IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        MESSAGE("-- MacOS")
        ADD_DEFINITIONS("-fast -msse3")
      ENDIF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  ENDIF (CMAKE_COMPILER_IS_GNUCXX)
ENDIF (CMAKE_BUILD_TYPE STREQUAL "Release")

#
# Packaging
#
SET(CPACK_GENERATOR ZIP)
SET(CPACK_PACKAGE_FILE_NAME "curve_extraction")
INCLUDE(CPack)

IF (${BUILD_EXAMPLES})
  MESSAGE("-- Build the examples.")
  ADD_SUBDIRECTORY(examples)
ELSE (${BUILD_EXAMPLES})
  MESSAGE("-- Do not build any example.")
ENDIF (${BUILD_EXAMPLES})

ADD_SUBDIRECTORY(data)
ADD_SUBDIRECTORY(source)
ADD_SUBDIRECTORY(tests)

