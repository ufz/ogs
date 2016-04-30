GET_SOURCE_FILES(SOURCES)
GET_SOURCE_FILES(SOURCES_IO IO)

set(SOURCES ${SOURCES} ${SOURCES_IO})

GET_SOURCE_FILES(SOURCES_IO_BASE_XML IO/XmlIO)
set(SOURCES ${SOURCES} ${SOURCES_IO_BASE_XML})

if(QT4_FOUND)
	GET_SOURCE_FILES(SOURCES_IO_QT_XML IO/XmlIO/Qt)
	set(SOURCES ${SOURCES} ${SOURCES_IO_QT_XML})
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

if(QT4_FOUND)
	target_link_libraries(BaseLib PUBLIC Qt4::QtXml Qt4::QtXmlPatterns)
	if(WIN32 AND CMAKE_CROSSCOMPILING AND OPENSSL_FOUND)
		target_link_libraries(BaseLib PUBLIC Qt4::QtNetwork ${OPENSSL_LIBRARIES} ws2_32)
	endif()
endif()

if(TARGET Eigen)
	add_dependencies(BaseLib Eigen)
endif()
if(TARGET Boost)
	add_dependencies(BaseLib Boost)
endif()
