INCLUDE(ExternalData)

SET(ExternalData_OBJECT_STORES "${ExternalData_OBJECT_STORES_DEFAULT}" CACHE STRING
	"Semicolon-separated list of local directories holding test data files in the layout %(algo)/%(hash).")
MARK_AS_ADVANCED(ExternalData_OBJECT_STORES)
IF(NOT ExternalData_OBJECT_STORES)
	SET(ExternalData_OBJECT_STORES "${CMAKE_SOURCE_DIR}/../ogs6-data")
	FILE(MAKE_DIRECTORY "${ExternalData_OBJECT_STORES}")
ENDIF()

SET(ExternalData_SOURCE_ROOT ${CMAKE_SOURCE_DIR}/Tests/Data)
SET(ExternalData_BINARY_ROOT ${CMAKE_BINARY_DIR}/Tests/Data)
SET(ExternalData_LINK_CONTENT MD5)

SET(ExternalData_URL_TEMPLATES "http://www.opengeosys.org/images/dev/%(algo)/%(hash)")

ADD_CUSTOM_TARGET(
	move-data
	COMMAND ${CMAKE_COMMAND}
	-DExternalData_SOURCE_ROOT=${ExternalData_SOURCE_ROOT}
	-DExternalData_LINK_CONTENT=${ExternalData_LINK_CONTENT}
	-DExternalData_OBJECT_STORES=${ExternalData_OBJECT_STORES}
	-P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/MoveDataToStore.cmake
	VERBATIM
)

ExternalData_Add_Test(
	data # data target specified above
	NAME External-Data-Test
	COMMAND ls -ll ${ExternalData_SOURCE_ROOT}/Elliptic/2d-quads-x1000-y1000/
	DATA{${ExternalData_SOURCE_ROOT}/Elliptic/2d-quads-x1000-y1000/,REGEX:.*}
)
