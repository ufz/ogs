ogs_add_library(ProcessLibBoundaryConditionPython
            PythonBoundaryCondition.cpp
            PythonBoundaryCondition.h
            PythonBoundaryConditionLocalAssembler.h
            PythonBoundaryConditionPythonSideInterface.h
            BHEInflowPythonBoundaryCondition.h
            BHEInflowPythonBoundaryConditionPythonSideInterface.h)

target_compile_definitions(ProcessLibBoundaryConditionPython PUBLIC
                           OGS_USE_PYTHON)

target_link_libraries(ProcessLibBoundaryConditionPython
                      PUBLIC BaseLib
                             MathLib
                             MeshLib
                             NumLib
                             ${Python3_LIBRARIES}
                      PRIVATE pybind11::pybind11)

# For the embedded Python module
ogs_add_library(ProcessLibBoundaryConditionPythonModule
            PythonBoundaryConditionModule.cpp PythonBoundaryConditionModule.h
            BHEInflowPythonBoundaryConditionModule.cpp
            BHEInflowPythonBoundaryConditionModule.h)

target_link_libraries(ProcessLibBoundaryConditionPythonModule
                      PUBLIC ProcessLibBoundaryConditionPython
                             pybind11::pybind11)
