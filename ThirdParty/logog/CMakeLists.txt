# logog needs thread support on linux
find_package( Threads )
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
target_include_directories(logog SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(BUILD_SHARED_LIBS)
    install(TARGETS logog LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
