if(IS_SUBPROJECT)
    include (CPack)
    return()
endif()

include(packaging/PackagingMacros)
include(packaging/ArchiveTestdata)

#### Packaging setup ####
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OGS-6 THM/C Simulator")
set(CPACK_PACKAGE_VENDOR "OpenGeoSys Community (http://www.opengeosys.org)")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "OGS-${OGS_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE.txt")
set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.md")
# set(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_SOURCE_DIR}/README.md")

# Package file name
if(OGS_BUILD_GUI)
    set(SUFFIX "-de")
endif()
if(OGS_BUILD_UTILS)
    set(SUFFIX "${SUFFIX}-utils")
endif()
if(OGS_USE_MPI)
    set(SUFFIX "${SUFFIX}-mpi")
endif()

if(APPLE)
    string(REGEX MATCH "(^[0-9]*)" TMP ${CMAKE_SYSTEM_VERSION})
    math(EXPR OSX_VERSION_MINOR "${CMAKE_MATCH_1} - 4")
    set(CPACK_PACKAGE_FILE_NAME
        "ogs-${OGS_VERSION}-OSX-10.${OSX_VERSION_MINOR}-x${BITS}${SUFFIX}")
    set(CPACK_SOURCE_PACKAGE_FILE_NAME ${CPACK_PACKAGE_FILE_NAME})
else()
    set(CPACK_PACKAGE_FILE_NAME "ogs-${OGS_VERSION}-${CMAKE_SYSTEM}-x${BITS}${SUFFIX}")
endif()

if (WIN32)
    include (packaging/PackagingWin)
endif()
if(UNIX)
    if(APPLE)
        include (packaging/PackagingMac)
    else()
        include (packaging/PackagingLinux)
    endif()
endif()

# Download additional content
if(OGS_DOWNLOAD_ADDITIONAL_CONTENT)
    DownloadAdditionalFilesForPackaging(
        URLS http://docs.opengeosys.org/assets/releases/head/docs/DataExplorer-Manual.pdf
        #     http://docs.opengeosys.org/assets/releases/head/docs/User-Manual.pdf
        #     http://docs.opengeosys.org/assets/releases/head/docs/Theory-Manual.pdf
        DESTINATION docs
        PACKAGE_GROUP ogs_docs
    )

    if(WIN32)
        DownloadAdditionalFilesForPackaging(
            URLS http://docs.opengeosys.org/assets/releases/head/win/OGSFileConverter.exe
            DESTINATION bin
            EXECUTABLE TRUE
            PACKAGE_GROUP ogs_converter
        )
    endif()
    if(APPLE)
        DownloadAdditionalFilesForPackaging(
            URLS http://docs.opengeosys.org/assets/releases/head/mac/OGSFileConverter
            DESTINATION bin
            EXECUTABLE TRUE
            PACKAGE_GROUP ogs_converter
        )
    endif()
endif()

include (CPack)

cpack_add_component_group(Applications
    DISPLAY_NAME Applications
    DESCRIPTION "OpenGeoSys applications"
    EXPANDED
    BOLD_TITLE
)

cpack_add_component_group(Utilities
    DISPLAY_NAME Utilities
    DESCRIPTION "OpenGeoSys utilities"
    EXPANDED
)

cpack_add_component(ogs_extras
    DISPLAY_NAME "Extra tools"
    DESCRIPTION "Miscellaneous tools."
    GROUP Utilities
)

cpack_add_component(ogs_docs
    DISPLAY_NAME "Documentation"
    DESCRIPTION "PDF documentation."
    GROUP Utilities
)


# Install shared libraries, copied to bin-dir from e.g Conan
foreach(PATTERN "*.dll" "*.dylib")
    file(GLOB MATCHED_FILES ${EXECUTABLE_OUTPUT_PATH}/${PATTERN})
    # message(STATUS ${MATCHED_FILES})
    install(FILES ${MATCHED_FILES} DESTINATION bin)
endforeach()

# macOS frameworks are directories, exclude header files
file(GLOB MATCHED_DIRECTORIES "${EXECUTABLE_OUTPUT_PATH}/*.framework")
install(DIRECTORY ${MATCHED_DIRECTORIES} DESTINATION bin
    PATTERN "Headers" EXCLUDE)
