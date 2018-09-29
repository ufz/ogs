if(APPLE)
    find_program(REALPATH_TOOL_PATH grealpath)
elseif(UNIX)
    find_program(REALPATH_TOOL_PATH realpath)
else()
    return()
endif()
find_program(ZIP_TOOL_PATH zip)
if(NOT REALPATH_TOOL_PATH OR NOT ZIP_TOOL_PATH)
    return()
endif()

add_custom_target(archive-data
    bash ${PROJECT_SOURCE_DIR}/scripts/packaging/archive-testdata.sh
    DEPENDS data
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Packaging testdata to ogs6-data.tar.gz and ogs6-data.zip" VERBATIM
)
