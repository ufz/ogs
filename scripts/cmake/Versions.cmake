include(${PROJECT_SOURCE_DIR}/ThirdParty/json-cmake/JSONParser.cmake)
file(READ ${PROJECT_SOURCE_DIR}/web/data/versions.json jsonFileString)
sbeParseJson(ogs jsonFileString)
# Provides variables, e.g. ogs.minimum_version.gcc
# Output all variables with
#   foreach(var ${ogs})
#     message("${var} = ${${var}}")
#   endforeach()
