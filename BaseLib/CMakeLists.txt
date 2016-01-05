GET_SOURCE_FILES(SOURCES)

list(APPEND SOURCES "${CMAKE_CURRENT_BINARY_DIR}/BuildInfo.cpp" BuildInfo.h)

# Create the library
add_library(BaseLib STATIC ${SOURCES})

set_target_properties(BaseLib PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(BaseLib INTERFACE
	logog
)

if(MSVC)
	target_link_libraries(BaseLib INTERFACE WinMM) # needed for timeGetTime
endif()

if(TARGET Eigen)
	add_dependencies(BaseLib Eigen)
endif()
if(TARGET Boost)
	add_dependencies(BaseLib Boost)
endif()
