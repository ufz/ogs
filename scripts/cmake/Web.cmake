if(NOT NPM OR NOT PIP OR NOT PYTHON_EXECUTABLE)
    return()
endif()

if(YARN)
    set(PACKAGE_MANAGER ${YARN})
else()
    set(PACKAGE_MANAGER "${NPM} install")
endif()

add_custom_target(web-install
    COMMAND ${PACKAGE_MANAGER}
    COMMAND ${PIP} install -r requirements.txt
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/web
    BYPRODUCTS ${CMAKE_SOURCE_DIR}/web/node_modules
)

add_custom_target(web-import
    COMMAND ${PYTHON_EXECUTABLE} import.py
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/web/import
    DEPENDS web-install
)

if(DEFINED OGS_WEB_BASE_URL)
    set(HUGO_ARGS --baseURL ${OGS_WEB_BASE_URL})
endif()

add_custom_target(web
    COMMAND ${NPM} run build -- ${HUGO_ARGS}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/web
    DEPENDS web-install web-import
)

add_custom_target(web-clean
    COMMAND ${NPM} run clean
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/web
)
