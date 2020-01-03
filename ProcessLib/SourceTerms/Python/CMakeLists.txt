add_library(ProcessLibSourceTermPython
            CreatePythonSourceTerm.cpp
            CreatePythonSourceTerm.h
            PythonSourceTerm.cpp
            PythonSourceTerm.h
            PythonSourceTermLocalAssembler.h
            PythonSourceTermPythonSideInterface.h)
if(BUILD_SHARED_LIBS)
    install(TARGETS ProcessLibSourceTermPython
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

target_compile_definitions(ProcessLibSourceTermPython PUBLIC OGS_USE_PYTHON)

target_link_libraries(ProcessLibSourceTermPython
                      PUBLIC BaseLib
                             MathLib
                             MeshLib
                             NumLib
                             logog
                             ${Python3_LIBRARIES}
                      PRIVATE pybind11::pybind11)

# For the embedded Python module
add_library(ProcessLibSourceTermPythonModule PythonSourceTermModule.cpp
                                             PythonSourceTermModule.h)
if(BUILD_SHARED_LIBS)
    install(TARGETS ProcessLibSourceTermPythonModule
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

target_link_libraries(ProcessLibSourceTermPythonModule
                      PUBLIC ProcessLibSourceTermPython pybind11::pybind11)
