# MKLFlags.cmake

function(SET_MKL_FLAGS TARGET)
    if (USE_MKL)
        message(STATUS "Configuring target ${TARGET} to use Intel MKL")

        target_compile_options(${TARGET} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
        target_include_directories(${TARGET} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
        target_link_libraries(${TARGET} PUBLIC $<LINK_ONLY:MKL::MKL>)
    else()
        message(STATUS "MKL not found or not enabled; ${TARGET} will not use Intel MKL")
    endif()
endfunction()