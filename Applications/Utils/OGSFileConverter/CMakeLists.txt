add_library(OGSFileConverterLib FileListDialog.h OGSFileConverter.h
                                FileListDialog.cpp OGSFileConverter.cpp)
target_link_libraries(OGSFileConverterLib
                      PUBLIC ApplicationsFileIO MathLib QtBase ${VTK_LIBRARIES}
                      INTERFACE MeshLib)
if(BUILD_SHARED_LIBS)
    install(TARGETS OGSFileConverterLib
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

set_target_properties(OGSFileConverterLib PROPERTIES AUTOMOC TRUE AUTOUIC TRUE)

add_executable(OGSFileConverter main.cpp)

target_link_libraries(OGSFileConverter
                      PUBLIC ApplicationsFileIO GitInfoLib OGSFileConverterLib)

add_autogen_include(OGSFileConverterLib)

set_target_properties(OGSFileConverter OGSFileConverterLib
                      PROPERTIES FOLDER "Utilities")

if(OGS_USE_PCH)
    cotire(OGSFileConverterLib)
endif()

install(TARGETS OGSFileConverter RUNTIME DESTINATION bin COMPONENT Utilities)
