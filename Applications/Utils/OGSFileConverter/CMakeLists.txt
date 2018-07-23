set(SOURCES
    FileListDialog.cpp
    OGSFileConverter.cpp
)

include_directories(
    # Qt generated file includes
    ${CMAKE_CURRENT_BINARY_DIR}

    # Workaround for CMake 3.8
    ${CMAKE_CURRENT_BINARY_DIR}/OGSFileConverterLib_autogen/include
)

add_library(OGSFileConverterLib
    ${SOURCES}
    ${HEADERS}
)
target_link_libraries(OGSFileConverterLib
    PUBLIC QtBase MathLib
    INTERFACE MeshLib ApplicationsFileIO
)

add_executable(OGSFileConverter main.cpp)

target_link_libraries(OGSFileConverter
    OGSFileConverterLib
    ApplicationsFileIO
)

set_target_properties(OGSFileConverter OGSFileConverterLib
    PROPERTIES
    FOLDER "Utilities"
    AUTOMOC ON
    AUTOUIC ON
)

if(OGS_USE_PCH)
    cotire(OGSFileConverter)
endif()

####################
### Installation ###
####################

install(TARGETS OGSFileConverter RUNTIME DESTINATION bin COMPONENT ogs_cli)
set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "OGSFileConverter" "File Converter")
set(CPACK_NSIS_MENU_LINKS ${CPACK_NSIS_MENU_LINKS} "bin/OGSFileConverter.exe" "File Converter" PARENT_SCOPE)
