add_library(NetCdfDialogLib NetCdfConfigureDialog.cpp NetCdfConfigureDialog.h)
if(BUILD_SHARED_LIBS)
    install(TARGETS NetCdfConfigureDialog LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(NetCdfDialogLib PUBLIC Qt5::Widgets PRIVATE MathLib vtknetcdfcpp)
set_property(TARGET NetCdfDialogLib PROPERTY FOLDER "DataExplorer")

# Workaround for system installed VTK (tested on arch)
if(NOT OGS_USE_CONAN)
    target_include_directories(NetCdfDialogLib SYSTEM
        PUBLIC ${VTK_INSTALL_PREFIX}/include/vtk/vtknetcdfcpp)
endif()
