add_library(NetCdfDialogLib NetCdfConfigureDialog.cpp NetCdfConfigureDialog.h)
target_link_libraries(NetCdfDialogLib PUBLIC Qt5::Widgets PRIVATE MathLib vtknetcdfcpp)
set_property(TARGET NetCdfDialogLib PROPERTY FOLDER "DataExplorer")
