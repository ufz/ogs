# Source files
GET_SOURCE_FILES(SOURCES)

# Library
add_library(InSituLib ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS InSituLib LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

target_link_libraries(PUBLIC BaseLib PRIVATE MeshLib)

if(OGS_USE_PCH)
    cotire(InSituLib)
endif()
