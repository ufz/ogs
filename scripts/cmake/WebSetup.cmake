set(OGS_WEB_WORKSPACE_DIR ${PROJECT_BINARY_DIR}/web)
set(OGS_WEB_GENERATED_CONTENT_MOUNTS
    "[[module.mounts]]\nsource = '${OGS_WEB_WORKSPACE_DIR}/content'\ntarget = 'content'\n")
set(OGS_WEB_GENERATED_DATA_MOUNTS
    "[[module.mounts]]\nsource = '${OGS_WEB_WORKSPACE_DIR}/data'\ntarget = 'data'\n")
set(OGS_WEB_DOXYGEN_URL "https://ogs.ogs.xyz/ogs/doxygen/")
if(DEFINED ENV{OGS_COMMIT_TAG})
    set(OGS_WEB_DOXYGEN_URL
        "https://doxygen.opengeosys.org/$ENV{OGS_COMMIT_TAG}/")
endif()

file(MAKE_DIRECTORY ${OGS_WEB_WORKSPACE_DIR} ${OGS_WEB_WORKSPACE_DIR}/content
                    ${OGS_WEB_WORKSPACE_DIR}/data ${OGS_WEB_WORKSPACE_DIR}/static)
configure_file(${PROJECT_SOURCE_DIR}/web/hugo.toml.in
               ${OGS_WEB_WORKSPACE_DIR}/hugo.toml @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/web/package.json.in
               ${OGS_WEB_WORKSPACE_DIR}/package.json @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/web/tailwind.config.js.in
               ${OGS_WEB_WORKSPACE_DIR}/tailwind.config.js @ONLY)

add_custom_target(web
    # Generates web/data/stats.toml (ctest count, ...) for the Hugo site. Run
    # via ctest -N instead of at configure time since tests registered by
    # add_subdirectory()s after include(WebSetup) wouldn't be counted yet.
    COMMAND ${CMAKE_COMMAND}
            -DCTEST_COMMAND=${CMAKE_CTEST_COMMAND}
            -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
            -DBINARY_DIR=${PROJECT_BINARY_DIR}
            -P ${PROJECT_SOURCE_DIR}/scripts/cmake/web/GenerateWebStatsData.cmake
    COMMAND yarn
    COMMAND yarn build
    WORKING_DIRECTORY ${OGS_WEB_WORKSPACE_DIR}
    USES_TERMINAL
    COMMENT "Building web site in ${OGS_WEB_WORKSPACE_DIR}/public"
)

add_custom_target(preview-web
    COMMAND yarn
    COMMAND yarn server
    WORKING_DIRECTORY ${OGS_WEB_WORKSPACE_DIR}
    USES_TERMINAL
)
add_dependencies(preview-web web)
