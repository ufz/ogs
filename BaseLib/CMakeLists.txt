GET_SOURCE_FILES(SOURCES)

list(APPEND SOURCES "${CMAKE_CURRENT_BINARY_DIR}/BuildInfo.cpp" BuildInfo.h)

# Create the library
add_library(BaseLib STATIC ${SOURCES})

set_target_properties(BaseLib PROPERTIES LINKER_LANGUAGE CXX)

include_directories(
	.
)

target_link_libraries(BaseLib
	${Boost_LIBRARIES}
)

if(WIN32)
	target_link_libraries(BaseLib WinMM)
endif()

if(TARGET Eigen)
	add_dependencies(BaseLib Eigen)
endif()
