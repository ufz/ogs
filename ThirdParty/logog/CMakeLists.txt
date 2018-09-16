cmake_minimum_required (VERSION 2.8.4)
enable_testing()
project (logog)
include(CTest)
# logog needs thread support on linux
find_package( Threads )
set( RUNTIME_OUTPUT_DIRECTORY bin/ )
set( ARCHIVE_OUTPUT_DIRECTORY bin/ )
set( LIBRARY_OUTPUT_DIRECTORY bin/ )
set( CMAKE_BINARY_DIR build/ )
set( CMAKE_LEGACY_CYGWIN_WIN32 0 )
if(MSVC)
  # Force to always compile with W4
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-variadic-macros -Wall -Wno-long-long")
endif()
include_directories( include )
add_library( logog
	src/api.cpp
	src/checkpoint.cpp
	src/formatter.cpp
	src/lobject.cpp
	src/lstring.cpp
	src/message.cpp
	src/mutex.cpp
	src/node.cpp
	src/platform.cpp
	src/socket.cpp
	src/statics.cpp
	src/target.cpp
	src/timer.cpp
	src/topic.cpp
	src/unittest.cpp
)
set_target_properties(logog PROPERTIES DEBUG_POSTFIX "d")
add_executable( test-logog test/test.cpp )
target_link_libraries( test-logog logog ${CMAKE_THREAD_LIBS_INIT})
add_test( NAME test-harness COMMAND test-logog )
if (CMAKE_PROJECT_NAME STREQUAL logog)
	install(TARGETS logog ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
	install(DIRECTORY include/ DESTINATION "${CMAKE_INSTALL_PREFIX}/include/logog"
          FILES_MATCHING PATTERN "*.hpp")

	# Docs generation with Doxygen
	find_package( Doxygen )
	if( DOXYGEN_FOUND )
		add_custom_target (logog-doc ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/doxyfile
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
			COMMENT "Generating source code documentation with Doxygen." VERBATIM)
    endif() # DOXYGEN_FOUND
else()
    if(BUILD_SHARED_LIBS)
        install(TARGETS logog LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
    endif()
endif () # CMAKE_PROJECT_NAME = logog
