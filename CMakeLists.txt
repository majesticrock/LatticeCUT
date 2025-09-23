cmake_minimum_required(VERSION 3.16)

project(LatticeCUT VERSION 2.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_PREFIX_PATH "~/usr/local")
include_directories(
    $ENV{HOME}/usr/local/include
)

set(MKL_ARCH intel64)
set(MKL_THREADING parallel)
set(MKL_INTERFACE lp64)
set(MKL_LINK static)

find_package(MKL CONFIG QUIET PATHS $ENV{MKLROOT})
if (MKL_FOUND)
    set(USE_MKL ON CACHE BOOL "Use Intel MKL" FORCE)
else()
    set(USE_MKL OFF CACHE BOOL "Use Intel MKL" FORCE)
endif()

include(${CMAKE_SOURCE_DIR}/cmake/MKLFlags.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/SelectClusterBuildOptions.cmake)

find_package(OpenMP REQUIRED)
find_package(mrock REQUIRED)

include(${CMAKE_CURRENT_SOURCE_DIR}/../../PhdUtility/cmake/mrock-information.cmake)

add_executable(latticecut 
    sources/LatticeCUT.cpp
)
target_sources(latticecut PRIVATE
    sources/LatticeCUT/DOSModel.cpp
    sources/LatticeCUT/ModeHelper.cpp
)
SET_COMPILER_FLAGS(latticecut)

add_subdirectory(sources/DOS)

target_link_libraries(latticecut PRIVATE OpenMP::OpenMP_CXX DOS ${mrock_LIBRARIES})
target_include_directories(latticecut PRIVATE ${mrock_INCLUDE_DIRS})