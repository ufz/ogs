GET_SOURCE_FILES(SOURCES)
APPEND_SOURCE_FILES(SOURCES IO)

APPEND_SOURCE_FILES(SOURCES IO/XmlIO)

if(Qt5XmlPatterns_FOUND)
    APPEND_SOURCE_FILES(SOURCES IO/XmlIO/Qt)
endif()

list(APPEND SOURCES "${CMAKE_CURRENT_BINARY_DIR}/BuildInfo.cpp" BuildInfo.h)

# Create the library
add_library(BaseLib ${SOURCES})

set_target_properties(BaseLib PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(BaseLib
    logog
)

if(MSVC)
    target_link_libraries(BaseLib WinMM) # needed for timeGetTime
endif()

if(Qt5XmlPatterns_FOUND)
    target_link_libraries(BaseLib Qt5::Xml Qt5::XmlPatterns)
    if(WIN32 AND CMAKE_CROSSCOMPILING AND OPENSSL_FOUND)
        target_link_libraries(BaseLib Qt5::Network ${OPENSSL_LIBRARIES} ws2_32)
    endif()
endif()

if(TARGET Eigen)
    add_dependencies(BaseLib Eigen)
endif()
if(TARGET Boost)
    add_dependencies(BaseLib Boost)
endif()

cotire(BaseLib)
