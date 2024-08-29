# cmake-lint: disable=C0103

# Prefer more recent Python version
set(Python_FIND_STRATEGY VERSION)

# Don't use venv
set(Python_FIND_VIRTUALENV STANDARD)

find_package(Python ${python_version} COMPONENTS Interpreter REQUIRED)

if(${Python_VERSION} VERSION_GREATER_EQUAL 3.9)
    set(_upgrade_deps --upgrade-deps)
endif()

if(UV_TOOL_PATH)
    set(_create_venv_command ${UV_TOOL_PATH} venv --python ${Python_EXECUTABLE})
    set(_venv_tool "uv")
else()
    set(_create_venv_command ${Python_EXECUTABLE} -m venv ${_upgrade_deps}
                             .venv
    )
    set(_venv_tool "pip")
endif()

message(STATUS "┌─ PythonCreateVirtualEnv.cmake (using ${_venv_tool})")
list(APPEND CMAKE_MESSAGE_INDENT "│    ")

execute_process(
    COMMAND ${_create_venv_command} WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
)

list(POP_BACK CMAKE_MESSAGE_INDENT)
message(STATUS "└─ End PythonCreateVirtualEnv.cmake")
