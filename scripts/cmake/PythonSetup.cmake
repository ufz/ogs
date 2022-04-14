# cmake-lint: disable=C0103

# Prefer unix location over frameworks (Apple-only)
set(Python_FIND_FRAMEWORK LAST)

# Prefer more recent Python version
set(Python_FIND_STRATEGY VERSION)
set(_python_componets Interpreter)
if(OGS_USE_PYTHON OR OGS_USE_PIP)
    set(CMAKE_REQUIRE_FIND_PACKAGE_Python TRUE)
endif()
if(OGS_USE_PYTHON)
    list(APPEND _python_componets Development)
endif()
find_package(
    Python ${ogs.minimum_version.python} COMPONENTS ${_python_componets}
)

if(OGS_USE_PIP)
    set(OGS_PYTHON_PACKAGES ""
        CACHE INTERNAL "List of Python packages to be installed via pip."
    )

    if(NOT EXISTS ${PROJECT_BINARY_DIR}/.venv)
        execute_process(
            COMMAND ${Python_EXECUTABLE} -m venv .venv
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
        unset(_OGS_PYTHON_PACKAGES_SHA1 CACHE)
    endif()
    set(Python_ROOT_DIR ${PROJECT_BINARY_DIR}/.venv)
    if(MSVC)
        set(Python_EXECUTABLE ${Python_ROOT_DIR}/Scripts/python.exe)
        set(LOCAL_VIRTUALENV_BIN_DIR ${PROJECT_BINARY_DIR}/.venv/Scripts
            CACHE INTERNAL ""
        )
    else()
        set(Python_EXECUTABLE ${Python_ROOT_DIR}/bin/python)
        set(LOCAL_VIRTUALENV_BIN_DIR ${PROJECT_BINARY_DIR}/.venv/bin
            CACHE INTERNAL ""
        )
    endif()
    if(OGS_BUILD_TESTING)
        # Notebook requirements from versions.json
        foreach(var ${ogs.python.notebook_requirements})
            list(APPEND OGS_PYTHON_PACKAGES
                 "${ogs.python.notebook_requirements_${var}}"
            )
        endforeach()
        if("${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}" VERSION_EQUAL
           "3.10"
        )
            # Default python on arch is 3.10; there are no vtk wheels for it
            # (VTUInterface). Use wheels for vtk from pyvista.
            list(APPEND OGS_PYTHON_PACKAGES
                 "--find-links https://wheels.pyvista.org vtk"
            )
        endif()
        list(APPEND OGS_PYTHON_PACKAGES
             "snakemake==${ogs.minimum_version.snakemake}"
        )
        set(SNAKEMAKE ${LOCAL_VIRTUALENV_BIN_DIR}/snakemake CACHE FILEPATH ""
                                                                  FORCE
        )
    endif()
endif()
