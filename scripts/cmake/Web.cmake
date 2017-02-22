if(NOT NPM)
    return()
endif()

if(YARN)
    set(PACKAGE_MANAGER ${YARN})
else()
    set(PACKAGE_MANAGER "${NPM} install")
endif()

add_custom_target(web-install
    COMMAND ${PACKAGE_MANAGER}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/web
    BYPRODUCTS ${CMAKE_SOURCE_DIR}/web/node_modules
)

if(PIP AND PYTHON_EXECUTABLE AND
    (EXISTS ${CMAKE_SOURCE_DIR}/web/import/secret.py OR
     DEFINED ENV{CONTENTFUL_ACCESS_TOKEN}))

    add_custom_target(web-import
        COMMAND ${PIP} install -r ../requirements.txt
        COMMAND ${PYTHON_EXECUTABLE} import.py
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/web/import
        DEPENDS web-install
    )

    set(IMPORT_TARGET web-import)
else()
    message(STATUS "[web] Skipping import from Contentful!")
endif()

if(DEFINED OGS_WEB_BASE_URL)
    set(HUGO_ARGS --baseURL ${OGS_WEB_BASE_URL})
endif()

if(DEFINED ENV{JENKINS_URL})
    set(HUGO_ARGS ${HUGO_ARGS} --canonifyURLs)
endif()

add_custom_target(web
    COMMAND ${NPM} run build:release -- ${HUGO_ARGS}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/web
    DEPENDS web-install ${IMPORT_TARGET}
)

add_custom_target(web-clean
    COMMAND ${NPM} run clean
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/web
)
