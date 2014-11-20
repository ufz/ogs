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

IF(HOSTNAME STREQUAL "envinf1.eve.ufz.de")
	ADD_CUSTOM_TARGET(
		sync-data
		COMMAND ${CMAKE_COMMAND} -E copy_directory
		${CMAKE_SOURCE_DIR}/../ogs6-data
		/data/ogs/ogs6-data
	)
	IF(CURL_TOOL_PATH)
		ADD_CUSTOM_COMMAND(
			TARGET sync-data POST_BUILD
			COMMAND ${CURL_TOOL_PATH} --insecure 'https://svn.ufz.de:8443/buildByToken/build?job=Tmp_Trigger&token=ogsbuild&cause=Triggered_by_sync-data_target_on_envinf1'
			COMMENT "Triggered sync to opengeosys.org, see https://svn.ufz.de:8443/job/OGS-6/job/SyncExternalData"
		)
	ELSE()
		MESSAGE(STATUS "curl tool was not found but is required for the sync-data target!")
	ENDIF()
ENDIF()
