if(BUILD_TESTING)
    CPMAddPackage(
        NAME googletest
        GITHUB_REPOSITORY google/googletest
        GIT_TAG 389cb68b87193358358ae87cc56d257fd0d80189
        OPTIONS
            "INSTALL_GTEST OFF"
            "gtest_force_shared_crt ON"
    )

    CPMAddPackage(
        NAME autocheck
        GITHUB_REPOSITORY ufz/autocheck
        GIT_TAG e388ecbb31c49fc2724c8d0436da313b6edca7fd
        DOWNLOAD_ONLY YES
    )
    if(autocheck_ADDED)
        add_library(autocheck INTERFACE IMPORTED)
        target_include_directories(autocheck INTERFACE ${autocheck_SOURCE_DIR}/include)
    endif()

    CPMAddPackage(
        NAME vtkdiff
        GITHUB_REPOSITORY ufz/vtkdiff
        GIT_TAG 49403cee266bb8e80405a02d677dbb5f71afc61a
        OPTIONS
            "VTK_LIBRARIES vtkIOXML"
    )
    if(vtkdiff_ADDED)
        install(PROGRAMS $<TARGET_FILE:vtkdiff> DESTINATION bin)
    endif()
endif()

CPMAddPackage(
    NAME exprtk
    GITHUB_REPOSITORY ArashPartow/exprtk
    GIT_TAG c7c219480d9678eec7383a4a99030683c4a84d91
    DOWNLOAD_ONLY YES
)
if(exprtk_ADDED)
    add_library(exprtk INTERFACE IMPORTED)
    target_include_directories(exprtk INTERFACE ${exprtk_SOURCE_DIR})
endif()

CPMAddPackage(
    NAME spdlog
    GITHUB_REPOSITORY gabime/spdlog
    VERSION 1.8.2
)

CPMAddPackage(
    NAME tclap
    GITHUB_REPOSITORY ufz/tclap
    GIT_TAG 03abc3a3327214137c6ffd5b9a6efe23f0927cc2
    DOWNLOAD_ONLY YES
)
if(tclap_ADDED)
    add_library(tclap INTERFACE IMPORTED)
    target_include_directories(tclap INTERFACE ${tclap_SOURCE_DIR}/include)
endif()

CPMAddPackage(
    NAME tetgen
    GITHUB_REPOSITORY ufz/tetgen
    GIT_TAG 603ba181ebfaed38eec88532720e282606009b73
)
if(tetgen_ADDED)
    install(PROGRAMS $<TARGET_FILE:tetgen> DESTINATION bin)
endif()

if(OGS_USE_PYTHON)
    CPMAddPackage(
        NAME pybind11
        GITHUB_REPOSITORY pybind/pybind11
        GIT_TAG f1abf5d9159b805674197f6bc443592e631c9130
        # pybind11 uses old CMake find functionality, pass variables to use
        # the same Python installation.
        OPTIONS
            "PYTHON_INCLUDE_DIR ${Python3_INCLUDE_DIRS}"
            "PYTHON_LIBRARIES ${Python3_LIBRARIES}"
            "PYTHON_EXECUTABLE ${Python3_EXECUTABLE}"
            "PYBIND11_PYTHON_VERSION ${Python3_VERSION}"
    )
endif()

if (OGS_BUILD_PROCESS_ComponentTransport
    OR OGS_BUILD_PROCESS_RichardsComponentTransport)
    CPMAddPackage(
        NAME iphreeqc
        GITHUB_REPOSITORY ufz/iphreeqc
        GIT_TAG b1047d3eb03e7ef1b850231be35acb9c6a2cf345
        DOWNLOAD_ONLY YES
    )
    if(iphreeqc_ADDED)
        include(scripts/cmake/iphreeqc.cmake)
    endif()
endif()
