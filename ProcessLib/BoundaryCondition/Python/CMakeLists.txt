ogs_add_library(
    ProcessLibBoundaryConditionPython
    PythonBoundaryCondition.cpp
    PythonBoundaryCondition.h
    PythonBoundaryConditionLocalAssembler.h
    PythonBoundaryConditionPythonSideInterface.h
    BHEInflowPythonBoundaryCondition.h
    BHEInflowPythonBoundaryConditionPythonSideInterface.h
)

target_compile_definitions(
    ProcessLibBoundaryConditionPython
    PUBLIC OGS_USE_PYTHON
           # Enabled elements
           OGS_MAX_ELEMENT_DIM=${OGS_MAX_ELEMENT_DIM}
           OGS_MAX_ELEMENT_ORDER=${OGS_MAX_ELEMENT_ORDER}
           $<$<BOOL:${OGS_ENABLE_ELEMENT_SIMPLEX}>:OGS_ENABLE_ELEMENT_SIMPLEX>
           $<$<BOOL:${OGS_ENABLE_ELEMENT_CUBOID}>:OGS_ENABLE_ELEMENT_CUBOID>
           $<$<BOOL:${OGS_ENABLE_ELEMENT_PRISM}>:OGS_ENABLE_ELEMENT_PRISM>
           $<$<BOOL:${OGS_ENABLE_ELEMENT_PYRAMID}>:OGS_ENABLE_ELEMENT_PYRAMID>
)

target_link_libraries(
    ProcessLibBoundaryConditionPython PUBLIC BaseLib MathLib MeshLib NumLib
                                             ${Python3_LIBRARIES}
    PRIVATE pybind11::pybind11
)

# For the embedded Python module
ogs_add_library(
    ProcessLibBoundaryConditionPythonModule PythonBoundaryConditionModule.cpp
    PythonBoundaryConditionModule.h BHEInflowPythonBoundaryConditionModule.cpp
    BHEInflowPythonBoundaryConditionModule.h
)

target_link_libraries(
    ProcessLibBoundaryConditionPythonModule
    PUBLIC ProcessLibBoundaryConditionPython pybind11::pybind11
)
