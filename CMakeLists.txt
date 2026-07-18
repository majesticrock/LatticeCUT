cmake_minimum_required(VERSION 3.30)

project(LatticeCUT VERSION 2.1.0 LANGUAGES CXX)

# Set local user-specific settings
# By default, this adds ~/usr/local/ and ../../.mrock/ to the cmake and include paths
# If additional configurations/includes are needed, provide them in cmake/user_settings.cmake
# Extra include dirs are stored in the cmake variable EXTRA_INCLUDE_DIRS
include(cmake/include_settings.cmake)

# Build type settings; defaults to Release
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "Build type is '${CMAKE_BUILD_TYPE}'.")

#
# Library specifics
#

find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(mrock 
    REQUIRED 
    COMPONENTS 
        symbolic_operators 
        iEoM 
        utility
)

add_executable(latticecut 
    sources/LatticeCUT.cpp
)
target_sources(latticecut PRIVATE
    sources/LatticeCUT/DOSModel.cpp
    sources/LatticeCUT/ModeHelper.cpp
    sources/LatticeCUT/T_C.cpp
    sources/LatticeCUT/ComplexModel.cpp
)

# provided by the mrock library
mrock_generate_information_header(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/generated/latticecut/info.h"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    OUT_HEADER MROCK_GENERATED_INFO_HEADER
    OUT_INCLUDE_DIR MROCK_GENERATED_INCLUDE_DIR
)
include(${CMAKE_CURRENT_LIST_DIR}/cmake/latticecut_options.cmake)
# Sets include dir for the generated metadata header file
target_include_directories(latticecut_options INTERFACE ${MROCK_GENERATED_INCLUDE_DIR})

option(LATTICE_CUT_RESIDUALS "Compute residuals and eigenvectors of the transformed dynamical matrix." OFF)
if (LATTICE_CUT_RESIDUALS)
    target_compile_definitions(latticecut PRIVATE LATTICE_CUT_RESIDUALS)
endif()
option(LATTICE_CUT_FULL_DIAG "Fully diagonalize the transformed dynamical matrix." OFF)
if (LATTICE_CUT_FULL_DIAG)
    target_compile_definitions(latticecut PRIVATE LATTICE_CUT_FULL_DIAG)
endif()

add_subdirectory(sources/DOS)

target_link_libraries(latticecut 
    PRIVATE 
        OpenMP::OpenMP_CXX 
        ${mrock_LIBRARIES}
        DOS
        latticecut_options
)

option(BUILD_TESTING 
    "Build and configure the testing suit" 
    OFF
)

if(BUILD_TESTING)
    include(tests/test_config.cmake)
endif()