# From http://www.cmake.org/pipermail/cmake/2012-September/052098.html
macro(ConfigureMacOSXBundle TARGET_NAME ICON_FILE_PATH)

    get_filename_component(ICON_FILE_NAME "${ICON_FILE_PATH}" NAME)

    set_target_properties(${TARGET_NAME} PROPERTIES
        MACOSX_BUNDLE_INFO_STRING "${PROJECT_NAME} ${OGS_VERSION} - ${TARGET_NAME}"
        MACOSX_BUNDLE_ICON_FILE ${ICON_FILE_NAME}
        MACOSX_BUNDLE_GUI_IDENTIFIER "org.opengeosys"
        MACOSX_BUNDLE_BUNDLE_NAME ${PROJECT_NAME}-${TARGET_NAME}
        MACOSX_BUNDLE_SHORT_VERSION_STRING ${OGS_VERSION}
        MACOSX_BUNDLE_LONG_VERSION_STRING "${PROJECT_NAME} ${OGS_VERSION}"
        MACOSX_BUNDLE_BUNDLE_VERSION ${OGS_VERSION}
        MACOSX_BUNDLE_COPYRIGHT "Copyright (c) 2012-2018, OpenGeoSys Community. All Rights Reserved."
    )

    set_source_files_properties(${ICON_FILE_PATH} PROPERTIES
        MACOSX_PACKAGE_LOCATION Resources)
endmacro()

#
# DownloadAdditionalFilesForPackaging
# -------
#
# Downloads files (into build/AdditionalContent) and packages them. Order of arguments can be arbitrary.
#
# AddTest(
#   URLS <multiple URLs>
#   DESTINATION <install directory>
#   PACKAGE_GROUP <name of the package to add these files to> # optional, defaults to ogs_extras
#   EXECUTABLE <TRUE or FALSE, are these files executables?>  # optional, defaults to FALSE
# )
function(DownloadAdditionalFilesForPackaging)

    # parse args
    set(options NONE)
    set(oneValueArgs DESTINATION EXECUTABLE PACKAGE_GROUP)
    set(multiValueArgs URLS)
    cmake_parse_arguments(DownloadAdditionalFilesForPackaging
        "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    # set defaults
    if(NOT DownloadAdditionalFilesForPackaging_EXECUTABLE)
        set(DownloadAdditionalFilesForPackaging_EXECUTABLE FALSE)
    endif()
    if(NOT DownloadAdditionalFilesForPackaging_PACKAGE_GROUP)
        set(DownloadAdditionalFilesForPackaging_PACKAGE_GROUP ogs_extras)
    endif()

    foreach(URL ${DownloadAdditionalFilesForPackaging_URLS})
        get_filename_component(FILE_NAME ${URL} NAME)
        get_filename_component(FILE_EXTENSION ${URL} EXT)
        set(FILE_PATH ${PROJECT_BINARY_DIR}/AdditionalContent/${FILE_NAME})
        if(NOT EXISTS ${FILE_PATH})
            file(DOWNLOAD ${URL} ${FILE_PATH} SHOW_PROGRESS)
        endif()
        if(EXE)
            install(PROGRAMS ${FILE_PATH} DESTINATION ${DownloadAdditionalFilesForPackaging_DESTINATION} COMPONENT ${DownloadAdditionalFilesForPackaging_PACKAGE_GROUP})
        else()
            install(FILES ${FILE_PATH} DESTINATION ${DownloadAdditionalFilesForPackaging_DESTINATION} COMPONENT ${DownloadAdditionalFilesForPackaging_PACKAGE_GROUP})
        endif()
    endforeach()

endfunction()

macro(install_qt5_plugin _qt_plugin_name _qt_plugins_var)
    get_target_property(_qt_plugin_path "${_qt_plugin_name}" LOCATION)
    if(EXISTS "${_qt_plugin_path}")
        get_filename_component(_qt_plugin_file "${_qt_plugin_path}" NAME)
        get_filename_component(_qt_plugin_type "${_qt_plugin_path}" PATH)
        get_filename_component(_qt_plugin_type "${_qt_plugin_type}" NAME)
        if(APPLE)
            set(_qt_plugin_dir "PlugIns")
        elseif(WIN32)
            set(_qt_plugin_dir "plugins")
        endif()
        set(_qt_plugin_dest "${_qt_plugin_dir}/${_qt_plugin_type}")
        install(FILES "${_qt_plugin_path}"
            DESTINATION "${_qt_plugin_dest}"
            ${COMPONENT})
        set(${_qt_plugins_var}
            "${${_qt_plugins_var}};\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${_qt_plugin_dest}/${_qt_plugin_file}")
    else()
        message(FATAL_ERROR "QT plugin ${_qt_plugin_name} not found")
    endif()
endmacro()

