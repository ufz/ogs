ogs_add_library(ProcessLibSourceTermPython
                CreatePythonSourceTerm.cpp
                CreatePythonSourceTerm.h
                PythonSourceTerm.cpp
                PythonSourceTerm.h
                PythonSourceTermLocalAssembler.h
                PythonSourceTermPythonSideInterface.h)

target_compile_definitions(ProcessLibSourceTermPython PUBLIC OGS_USE_PYTHON)

target_link_libraries(ProcessLibSourceTermPython
                      PUBLIC BaseLib
                             MathLib
                             MeshLib
                             NumLib
                             ${Python3_LIBRARIES}
                      PRIVATE pybind11::pybind11)

# For the embedded Python module
ogs_add_library(ProcessLibSourceTermPythonModule PythonSourceTermModule.cpp
                                                 PythonSourceTermModule.h)

target_link_libraries(ProcessLibSourceTermPythonModule
                      PUBLIC ProcessLibSourceTermPython pybind11::pybind11)
