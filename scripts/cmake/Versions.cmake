CPMAddPackage(
    NAME json-cmake
    GITHUB_REPOSITORY ufz/json-cmake
    GIT_TAG 9708cb091f6b89b94d71ae98f8b9e68ea04c47dd
    DOWNLOAD_ONLY YES
)
include("${json-cmake_SOURCE_DIR}/JSONParser.cmake")
file(READ ${PROJECT_SOURCE_DIR}/web/data/versions.json jsonFileString)
sbeParseJson(ogs jsonFileString)
# Provides variables, e.g. ogs.minimum_version.gcc
# Output all variables with
#   foreach(var ${ogs})
#     message("${var} = ${${var}}")
#   endforeach()
