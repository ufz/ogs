ogs_add_library(
    ProcessLibSourceTermPython
    CreatePythonSourceTerm.cpp
    CreatePythonSourceTerm.h
    PythonSourceTerm.cpp
    PythonSourceTerm.h
    PythonSourceTermLocalAssembler.h
    PythonSourceTermPythonSideInterface.h
)

target_compile_definitions(
    ProcessLibSourceTermPython
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
    ProcessLibSourceTermPython PUBLIC BaseLib MathLib MeshLib NumLib
                                      ${Python3_LIBRARIES}
    PRIVATE pybind11::pybind11 spdlog::spdlog
)

# For the embedded Python module
ogs_add_library(
    ProcessLibSourceTermPythonModule PythonSourceTermModule.cpp
    PythonSourceTermModule.h
)

target_link_libraries(
    ProcessLibSourceTermPythonModule PUBLIC ProcessLibSourceTermPython
                                            pybind11::pybind11
)
