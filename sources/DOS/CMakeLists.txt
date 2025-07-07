set(DOS_SOURCES
    Base.cpp
    Selector.cpp
    FreeElectrons.cpp
    SimpleCubic.cpp
)

add_library(DOS STATIC ${DOS_SOURCES})

SET_COMPILER_FLAGS(DOS)

target_include_directories(DOS PUBLIC Boost::headers ${CMAKE_CURRENT_LIST_DIR})
