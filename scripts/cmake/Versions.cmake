CPMAddPackage(
    NAME json-cmake
    GITHUB_REPOSITORY ufz/json-cmake
    GIT_TAG 9708cb091f6b89b94d71ae98f8b9e68ea04c47dd
    DOWNLOAD_ONLY YES
)
include("${json-cmake_SOURCE_DIR}/JSONParser.cmake")
file(READ ${PROJECT_SOURCE_DIR}/web/data/versions.json jsonFileString)
set_property(
    DIRECTORY . APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
                                "${PROJECT_SOURCE_DIR}/web/data/versions.json"
)
sbeParseJson(ogs jsonFileString)
# Provides variables, e.g. ogs.minimum_version.gcc. Output all variables with
# ~~~
#   foreach(var ${ogs})
#     message("${var} = ${${var}}")
#   endforeach()
# ~~~

# Overwriting individual entries
set(OGS_DEPENDENCY_VERSIONS ""
    CACHE STRING "Overwrite versions.json; syntax: 'key=value;key=value'"
)

foreach(version ${OGS_DEPENDENCY_VERSIONS})
    if("${version}" MATCHES "^(.*)=(.*)$")
        message(
            STATUS
                "Overwriting versions.json: ogs.${CMAKE_MATCH_1} = ${CMAKE_MATCH_2} (default ${ogs.${CMAKE_MATCH_1}})"
        )
        # cmake-lint: disable=C0103
        set(ogs.${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
    endif()
endforeach()
