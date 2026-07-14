set(DOS_SOURCES
    Base.cpp
    Selector.cpp
    FreeElectrons.cpp
    SimpleCubic.cpp
    FCC.cpp
    BCC.cpp
    HoneyComb.cpp
    SinglePeak.cpp
    DoublePeak.cpp
    EnergyRanges.cpp
)

add_library(DOS STATIC ${DOS_SOURCES})

SET_COMPILER_FLAGS(DOS)

target_include_directories(DOS PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(DOS PRIVATE mrock::utility)