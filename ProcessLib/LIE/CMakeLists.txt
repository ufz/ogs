add_subdirectory(Common)
get_source_files(SOURCES HydroMechanics)
append_source_files(SOURCES HydroMechanics/LocalAssembler)
append_source_files(SOURCES SmallDeformation)
append_source_files(SOURCES SmallDeformation/LocalAssembler)

ogs_add_library(LIE ${SOURCES})
target_link_libraries(LIE PUBLIC ProcessLib LIECommon PRIVATE ParameterLib)

if(BUILD_TESTING)
    include(HydroMechanics/Tests.cmake)
    include(SmallDeformation/Tests.cmake)
endif()
