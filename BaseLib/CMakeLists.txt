GET_SOURCE_FILES(SOURCES)
GET_SOURCE_FILES(SOURCES_IO IO)

set(SOURCES ${SOURCES} ${SOURCES_IO})

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

if(TARGET Eigen)
	add_dependencies(BaseLib Eigen)
endif()
if(TARGET Boost)
	add_dependencies(BaseLib Boost)
endif()
