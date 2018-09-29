# Source files
GET_SOURCE_FILES(SOURCES_APPLICATIONSLIB)
set(LIB_SOURCES ${SOURCES_APPLICATIONSLIB})

# Library
add_library(ApplicationsLib ${LIB_SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS ApplicationsLib LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

target_link_libraries(ApplicationsLib
    PUBLIC BaseLib GeoLib Processes logog
    PRIVATE MathLib MeshLib MeshGeoToolsLib NumLib
)

# Set cpp definitions if the cmake option is enabled for the given process.
foreach(process ${ProcessesList})
    if(OGS_BUILD_PROCESS_${process})
        string(TOUPPER "OGS_BUILD_PROCESS_${process}" EnableProcess)
        set_property(SOURCE ProjectData.cpp APPEND PROPERTY COMPILE_DEFINITIONS
            ${EnableProcess})
    endif()
endforeach()

if(OGS_USE_PCH)
    cotire(ApplicationsLib)
endif()

if(OGS_USE_PYTHON)
    target_link_libraries(ApplicationsLib PRIVATE pybind11::pybind11)
endif()
