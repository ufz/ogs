add_library(ProcessLibBoundaryConditionPython
    PythonBoundaryCondition.cpp
    PythonBoundaryCondition.h
    PythonBoundaryConditionLocalAssembler.h
    PythonBoundaryConditionPythonSideInterface.h)

target_compile_definitions(ProcessLibBoundaryConditionPython
    PUBLIC OGS_USE_PYTHON)

target_link_libraries(ProcessLibBoundaryConditionPython
    PUBLIC BaseLib MathLib MeshLib NumLib logog
    PRIVATE pybind11::pybind11)

# For the embedded Python module
add_library(ProcessLibBoundaryConditionPythonModule
    PythonBoundaryConditionModule.cpp
    PythonBoundaryConditionModule.h)

target_link_libraries(ProcessLibBoundaryConditionPythonModule
    PUBLIC
    ProcessLibBoundaryConditionPython
    pybind11::pybind11)
