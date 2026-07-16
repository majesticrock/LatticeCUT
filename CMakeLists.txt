cmake_minimum_required(VERSION 3.30)

project(LatticeCUT VERSION 2.1.0 LANGUAGES CXX)

list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/usr/local")

find_package(OpenMP REQUIRED)
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