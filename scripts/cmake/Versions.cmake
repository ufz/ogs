include(${PROJECT_SOURCE_DIR}/scripts/cmake/json-cmake/JSONParser.cmake)
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
