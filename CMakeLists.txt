cmake_minimum_required(VERSION 3.0)
list(APPEND CMAKE_MESSAGE_CONTEXT star)
project(star VERSION 0.1 LANGUAGES CXX C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Include useful modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(CMakePrintHelpers)

# Add CPM Dependency Manager
include(FindCPM)

#############################################
# Options
#############################################

set(STAR_DEFAULT_BUILD_TYPE "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type specified. Setting CMAKE_BUILD_TYPE to ${STAR_DEFAULT_BUILD_TYPE}")
  set(CMAKE_BUILD_TYPE ${STAR_DEFAULT_BUILD_TYPE} CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

# Enable testing
option(STAR_BUILD_TESTS "Build tests for star" ON)

# Code Coverage
option(SLAP_CODE_COVERAGE "Compile star with Code Coverage." OFF)

# Documentation
option(STAR_BUILD_DOCS "Build documentation for star" OFF)

# Floating point precsion
set(STAR_FLOAT double CACHE STRING "Floating point precision for star (float, double)")

# Build with -march=native
option(STAR_VECTORIZE "Compile with -march=native" OFF)

##############################
# Dependencies
##############################
if (STAR_BUILD_DOCS)
  find_package(Doxygen REQUIRED)
endif ()
if (STAR_BUILD_TESTS)
  # Google Test (Testing)
  CPMAddPackage(
    NAME googletest
    GITHUB_REPOSITORY google/googletest
    VERSION 1.13.0
    DOWNLOAD_EXTRACT_TIMESTAMP
    OPTIONS
    "INSTALL_GTEST OFF"
    "gtest_force_shared_crt ON"
  )
  add_library(gtest::gtest ALIAS gtest_main)
  enable_testing()
  include(GoogleTest)
  include(CTest)
endif ()

#############################################
# Build
#############################################

# Add compile options
if (NOT WIN32)
  message(STATUS "Adding compile options")
  add_compile_options(-Wall -Wextra -pedantic -Werror -save-temps )
  if (STAR_VECTORIZE)
    add_compile_options(-march=native)
  endif ()
endif ()

add_subdirectory(src/star)

#############################################
# TESTING
#############################################
if (STAR_BUILD_TESTS)
  add_subdirectory(test)
endif()
