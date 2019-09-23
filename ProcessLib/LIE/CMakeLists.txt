add_subdirectory(Common)
append_source_files(SOURCES HydroMechanics)
append_source_files(SOURCES HydroMechanics/LocalAssembler)
append_source_files(SOURCES SmallDeformation)
append_source_files(SOURCES SmallDeformation/LocalAssembler)

add_library(LIE ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS LIE LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(LIE PUBLIC ProcessLib LIECommon PRIVATE ParameterLib)

include(HydroMechanics/Tests.cmake)
include(SmallDeformation/Tests.cmake)
