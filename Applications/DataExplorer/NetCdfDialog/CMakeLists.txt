add_library(NetCdfDialogLib NetCdfConfigureDialog.cpp NetCdfConfigureDialog.h)
if(BUILD_SHARED_LIBS)
    install(TARGETS NetCdfConfigureDialog
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(NetCdfDialogLib
                      PUBLIC Qt5::Widgets
                             ${NETCDF_LIBRARIES_CXX}
                             ${NETCDF_LIBRARIES_C}
                             ${HDF5_HL_LIBRARIES}
                             ${HDF5_LIBRARIES}
                      PRIVATE MathLib)
set_property(TARGET NetCdfDialogLib PROPERTY FOLDER "DataExplorer")

# Workaround for system installed VTK (tested on arch)
if(NOT OGS_USE_CONAN)
    target_include_directories(NetCdfDialogLib SYSTEM
                               PUBLIC ${NETCDF_INCLUDES_C}
                               PUBLIC ${NETCDF_INCLUDES_CXX})
endif()

add_autogen_include(NetCdfDialogLib)
