include(ExternalData)

set(ExternalData_LOCAL_FOLDER ogs6-data)

set(ExternalData_OBJECT_STORES "${ExternalData_OBJECT_STORES_DEFAULT}" CACHE STRING
	"Semicolon-separated list of local directories holding test data files in the layout %(algo)/%(hash).")
mark_as_advanced(ExternalData_OBJECT_STORES)
if(NOT ExternalData_OBJECT_STORES)
	set(ExternalData_OBJECT_STORES "${CMAKE_SOURCE_DIR}/../${ExternalData_LOCAL_FOLDER}")
	file(MAKE_DIRECTORY "${ExternalData_OBJECT_STORES}")
endif()

set(ExternalData_SOURCE_ROOT ${CMAKE_SOURCE_DIR}/Tests/Data)
set(ExternalData_BINARY_ROOT ${CMAKE_BINARY_DIR}/Tests/Data)
set(ExternalData_LINK_CONTENT MD5)

# Amazon S3 config
set(ExternalData_S3_Bucket opengeosys)

set(ExternalData_URL_TEMPLATES
	"http://www.opengeosys.org/images/dev/%(algo)/%(hash)"
	"http://${ExternalData_S3_Bucket}.s3.amazonaws.com/${ExternalData_LOCAL_FOLDER}/%(algo)/%(hash)"
)

add_custom_target(
	move-data
	COMMAND ${CMAKE_COMMAND}
	-DExternalData_SOURCE_ROOT=${ExternalData_SOURCE_ROOT}
	-DExternalData_LINK_CONTENT=${ExternalData_LINK_CONTENT}
	-DExternalData_OBJECT_STORES=${ExternalData_OBJECT_STORES}
	-P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/MoveDataToStore.cmake
	VERBATIM
)

if(S3CMD_TOOL_PATH)
	add_custom_target(
		upload-data
		COMMAND ${S3CMD_TOOL_PATH} sync --acl-public --skip-existing --no-check-md5
			${CMAKE_SOURCE_DIR}/../${ExternalData_LOCAL_FOLDER} s3://${ExternalData_S3_Bucket}
	)
endif()

if(HOSTNAME STREQUAL "envinf1.eve.ufz.de")
	add_custom_target(
		sync-data
		COMMAND ${CMAKE_COMMAND} -E copy_directory
		${CMAKE_SOURCE_DIR}/../${ExternalData_LOCAL_FOLDER}
		/data/ogs/${ExternalData_LOCAL_FOLDER}
	)
	if(CURL_TOOL_PATH)
		add_custom_command(
			TARGET sync-data POST_BUILD
			COMMAND ${CURL_TOOL_PATH} --insecure 'https://svn.ufz.de:8443/buildByToken/build?job=Tmp_Trigger&token=ogsbuild&cause=Triggered_by_sync-data_target_on_envinf1'
			COMMENT "Triggered sync to opengeosys.org, see https://svn.ufz.de:8443/job/OGS-6/job/SyncExternalData"
		)
	else()
		message(STATUS "curl tool was not found but is required for the sync-data target!")
	endif()
endif()
