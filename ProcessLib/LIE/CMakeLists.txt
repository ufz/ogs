APPEND_SOURCE_FILES(SOURCES BoundaryCondition)
APPEND_SOURCE_FILES(SOURCES Common)
APPEND_SOURCE_FILES(SOURCES HydroMechanics)
APPEND_SOURCE_FILES(SOURCES HydroMechanics/LocalAssembler)
APPEND_SOURCE_FILES(SOURCES SmallDeformation)
APPEND_SOURCE_FILES(SOURCES SmallDeformation/LocalAssembler)

add_library(LIE ${SOURCES})
target_link_libraries(LIE PUBLIC ProcessLib)

include(HydroMechanics/Tests.cmake)
include(SmallDeformation/Tests.cmake)
