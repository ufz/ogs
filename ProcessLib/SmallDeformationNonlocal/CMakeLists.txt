APPEND_SOURCE_FILES(SOURCES)

add_library(SmallDeformationNonlocal ${SOURCES})
target_link_libraries(SmallDeformationNonlocal
    PUBLIC ProcessLib
    PRIVATE ParameterLib
)
if(BUILD_SHARED_LIBS)
    install(TARGETS SmallDeformationNonlocal LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

include(Tests.cmake)
