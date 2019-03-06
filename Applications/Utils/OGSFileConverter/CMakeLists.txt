add_library(OGSFileConverterLib
    FileListDialog.h
    OGSFileConverter.h
    FileListDialog.cpp
    OGSFileConverter.cpp
)
target_link_libraries(OGSFileConverterLib
    PUBLIC ApplicationsFileIO MathLib QtBase ${VTK_LIBRARIES}
    INTERFACE MeshLib
)
if(BUILD_SHARED_LIBS)
    install(TARGETS OGSFileConverterLib LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

target_include_directories(OGSFileConverterLib PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/OGSFileConverterLib_autogen/include)

set_target_properties(OGSFileConverterLib PROPERTIES
    AUTOMOC TRUE
    AUTOUIC TRUE
)

add_executable(OGSFileConverter main.cpp)

target_link_libraries(OGSFileConverter
    PUBLIC OGSFileConverterLib ApplicationsFileIO
)

set_target_properties(OGSFileConverter OGSFileConverterLib
    PROPERTIES FOLDER "Utilities"
)

if(OGS_USE_PCH)
    cotire(OGSFileConverterLib)
endif()

install(TARGETS OGSFileConverter RUNTIME DESTINATION bin COMPONENT Utilities)
