ogs_add_library(
    ProcessLibBoundaryConditionAndSourceTermPython
    PythonBoundaryCondition.cpp
    PythonBoundaryCondition.h
    PythonBoundaryConditionLocalAssemblerInterface.h
    PythonBoundaryConditionLocalAssembler.h
    PythonBoundaryConditionPythonSideInterface.h
    BHEInflowPythonBoundaryCondition.h
    BHEInflowPythonBoundaryConditionPythonSideInterface.h
    CreatePythonSourceTerm.cpp
    CreatePythonSourceTerm.h
    PythonSourceTerm.cpp
    PythonSourceTerm.h
    PythonSourceTermLocalAssembler.h
    PythonSourceTermPythonSideInterface.h
    Utils/BcAndStLocalAssemblerImpl.h
    Utils/BcOrStData.h
    Utils/CollectAndInterpolateNodalDof.cpp
    Utils/CollectAndInterpolateNodalDof.h
    Utils/CreateLocalAssemblers.h
    Utils/LocalAssemblerFactoryPython.h
    Utils/NsAndWeight.h
)

target_compile_definitions(
    ProcessLibBoundaryConditionAndSourceTermPython
    PUBLIC # Enabled elements
           OGS_MAX_ELEMENT_DIM=${OGS_MAX_ELEMENT_DIM}
           OGS_MAX_ELEMENT_ORDER=${OGS_MAX_ELEMENT_ORDER}
           $<$<BOOL:${OGS_ENABLE_ELEMENT_SIMPLEX}>:OGS_ENABLE_ELEMENT_SIMPLEX>
           $<$<BOOL:${OGS_ENABLE_ELEMENT_CUBOID}>:OGS_ENABLE_ELEMENT_CUBOID>
           $<$<BOOL:${OGS_ENABLE_ELEMENT_PRISM}>:OGS_ENABLE_ELEMENT_PRISM>
           $<$<BOOL:${OGS_ENABLE_ELEMENT_PYRAMID}>:OGS_ENABLE_ELEMENT_PYRAMID>
)

target_link_libraries(
    ProcessLibBoundaryConditionAndSourceTermPython
    PUBLIC BaseLib MathLib MeshLib NumLib pybind11::module
)

if(BUILD_SHARED_LIBS AND MSVC)
    target_link_options(
        ProcessLibBoundaryConditionAndSourceTermPython PRIVATE
        /FORCE:UNRESOLVED
    )
endif()

# For the embedded Python module
ogs_add_library(
    ProcessLibBoundaryConditionAndSourceTermPythonModule
    PythonBoundaryConditionModule.cpp
    PythonBoundaryConditionModule.h
    BHEInflowPythonBoundaryConditionModule.cpp
    BHEInflowPythonBoundaryConditionModule.h
    PythonSourceTermModule.cpp
    PythonSourceTermModule.h
)

target_link_libraries(
    ProcessLibBoundaryConditionAndSourceTermPythonModule
    PUBLIC ProcessLibBoundaryConditionAndSourceTermPython
)
