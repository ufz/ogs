
add_executable(ogs ogs.cpp)

target_link_libraries(ogs
    PRIVATE BaseLib ApplicationsLib
)

if(OGS_USE_PYTHON)
    # Troubleshooting:
    # If you get linker errors, such as   ogs.cpp:(.text+0xb4): Warnung: undefinierter Verweis auf »_Py_ZeroStruct«
    # it could be that OGS is compiled with the wrong Python version.
    # I (Ch. Leh.) observed the following: The symbol _Py_ZeroStruct could not be found in /usr/lib/libpython3.6m.so (I intended to compile OGS with Python3).
    # It's apparently a Python2 symbol (present in /usr/lib/libpython2.7.so)
    # The compiler command-line was the following:
    # /usr/bin/g++ ... -DvtkRenderingVolume_AUTOINIT="1(vtkRenderingVolumeOpenGL2)" \
    #     -I/usr/include/vtk -I/usr/include/python2.7 -I/usr/include/freetype2 \
    #     -I/usr/include/libxml2 ... -I/.../BaseLib ... \
    #     -isystem /usr/include/python3.6m ... -o CMakeFiles/ogs.dir/ogs.cpp.o \
    #     -c /.../Applications/CLI/ogs.cpp
    # In particular, the Python2 include path comes before the Python3 include path.
    # Compiling OGS with Python2 solved the issue.
    # I assume (this is only a guess!) that VTK pulls in Python2 dependencies (on my system).
    # I assume that this is related to https://github.com/ufz/ogs/pull/2158.
    # Workaround: Always make sure that OGS is compiled with the same Python version as VTK.
    # TODO: Find out how to properly address the issue.

    # Performance warning from
    # https://github.com/pybind/pybind11/blob/master/docs/compiling.rst:
    # Since pybind11 is a metatemplate library, it is crucial that certain compiler
    # flags are provided to ensure high quality code generation. In contrast to the
    # pybind11_add_module() command, the CMake interface library only provides the
    # minimal set of parameters to ensure that the code using pybind11 compiles, but
    # it does not pass these extra compiler flags (i.e. this is up to you).
    # TODO: Enable further compiler/linker flags.
    target_link_libraries(ogs PRIVATE pybind11::embed)

    add_library(ogs_python_bindings STATIC ogs_python_bindings.cpp)
    target_link_libraries(ogs_python_bindings PRIVATE pybind11::embed)
    target_include_directories(ogs_python_bindings PRIVATE ${PYTHON_INCLUDE_DIRS})

    target_link_libraries(ogs PRIVATE ogs_python_bindings)

    if(BUILD_SHARED_LIBS)
        target_compile_definitions(ogs_python_bindings OGS_BUILD_SHARED_LIBS)
    endif()
endif()

if(OGS_USE_PETSC)
    target_link_libraries(ogs PRIVATE ${PETSC_LIBRARIES})
endif()

if(OGS_INSITU)
    target_link_libraries(ogs PRIVATE InSituLib)
endif()

if(OGS_USE_PCH)
    cotire(ogs)
endif()

####################
### Tests ##########
####################
add_test(NAME ogs_no_args COMMAND ogs)
set_tests_properties(ogs_no_args PROPERTIES WILL_FAIL TRUE)

####################
### Installation ###
####################
install(TARGETS ogs RUNTIME DESTINATION bin COMPONENT ogs_cli)

set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "ogs" "OGS Simulator")
cpack_add_component(ogs_cli
    DISPLAY_NAME "OGS THMC Simulator"
    DESCRIPTION "The command line interface for OpenGeoSys."
    GROUP Applications
)
