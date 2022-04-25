if(OGS_USE_PYTHON)
    # Troubleshooting: If you get linker errors, such as   ogs.cpp:(.text+0xb4):
    # undefined reference to `_Py_ZeroStruct' it could be that OGS is compiled
    # with the wrong Python version. I (Ch. Leh.) observed the following: The
    # symbol _Py_ZeroStruct could not be found in /usr/lib/libpython3.6m.so (I
    # intended to compile OGS with Python3). It's apparently a Python2 symbol
    # (present in /usr/lib/libpython2.7.so) The compiler command-line was the
    # following:
    # ~~~
    # /usr/bin/g++ ... -DvtkRenderingVolume_AUTOINIT="1(vtkRenderingVolumeOpenGL2)" \
    #     -I/usr/include/vtk -I/usr/include/python2.7 -I/usr/include/freetype2 \
    #     -I/usr/include/libxml2 ... -I/.../BaseLib ... \
    #     -isystem /usr/include/python3.6m ... -o CMakeFiles/ogs.dir/ogs.cpp.o \
    #     -c /.../Applications/CLI/ogs.cpp
    # ~~~
    # In particular, the Python2 include path comes before the Python3 include
    # path. Compiling OGS with Python2 solved the issue. I assume (this is only
    # a guess!) that VTK pulls in Python2 dependencies (on my system). I assume
    # that this is related to https://github.com/ufz/ogs/pull/2158. Workaround:
    # Always make sure that OGS is compiled with the same Python version as VTK.
    # The error described above should be detected automatically by cmake and an
    # appropriate message should be presented. The note is kept for the case
    # that the automatic detection does not work due to whatever reason.

    ogs_add_library(ogs_embedded_python ogs_embedded_python.cpp)

    # Performance warning from
    # https://github.com/pybind/pybind11/blob/master/docs/compiling.rst: Since
    # pybind11 is a metatemplate library, it is crucial that certain compiler
    # flags are provided to ensure high quality code generation. In contrast to
    # the pybind11_add_module() command, the CMake interface library only
    # provides the minimal set of parameters to ensure that the code using
    # pybind11 compiles, but it does not pass these extra compiler flags (i.e.
    # this is up to you). TODO: Enable further compiler/linker flags.
    target_link_libraries(
        ogs_embedded_python PUBLIC pybind11::embed
        PRIVATE ProcessLibBoundaryConditionAndSourceTermPythonModule
    )
    target_compile_definitions(
        ogs_embedded_python
        PUBLIC OGS_USE_PYTHON
               # Add macro definition, because static libs make special handling
               # necessary s.t. the embedded OpenGeoSys Python module won't be
               # removed by the linker.
               $<$<BOOL:${BUILD_SHARED_LIBS}>:PRIVATE
               OGS_BUILD_SHARED_LIBS>
    )

    # Create OpenGeoSys python module
    # https://pybind11.readthedocs.io/en/stable/compiling.html#building-with-cmake
    if(OGS_BUILD_PYTHON_MODULE)
        pybind11_add_module(
            OpenGeoSys MODULE ogs_python_module.cpp
            CommandLineArgumentParser.cpp
        )

        # lld linker strips out PyInit_OpenGeoSys symbol. Use standard linker.
        get_target_property(_link_options OpenGeoSys LINK_OPTIONS)
        if(_link_options)
            list(REMOVE_ITEM _link_options -fuse-ld=lld)
            set_target_properties(
                OpenGeoSys PROPERTIES LINK_OPTIONS "${_link_options}"
            )
        endif()

        target_link_libraries(
            OpenGeoSys PRIVATE ApplicationsLib BaseLib GitInfoLib tclap
                               pybind11::pybind11
        )
    endif()
endif()

ogs_add_executable(ogs ogs.cpp CommandLineArgumentParser.cpp)

target_link_libraries(
    ogs
    PRIVATE ApplicationsLib
            BaseLib
            CMakeInfoLib
            GitInfoLib
            MeshLib
            ProcessLib
            $<$<TARGET_EXISTS:ogs_embedded_python>:ogs_embedded_python>
            $<$<TARGET_EXISTS:MPI::MPI_CXX>:MPI::MPI_CXX>
            $<$<TARGET_EXISTS:InSituLib>:InSituLib>
            $<$<TARGET_EXISTS:petsc>:petsc>
            $<$<TARGET_EXISTS:VTK::ParallelMPI>:VTK::ParallelMPI>
            tclap
)

target_compile_definitions(ogs PRIVATE $<$<BOOL:${USE_INSITU}>:USE_INSITU>)

# ---- Tests ----
add_test(NAME ogs_no_args COMMAND ogs)
set_tests_properties(ogs_no_args PROPERTIES WILL_FAIL TRUE LABELS "default")

if(OGS_BUILD_PYTHON_MODULE AND NOT (WIN32 AND "${CMAKE_BUILD_TYPE}" STREQUAL
                                              "Release")
)
    # TODO: does not work on Windows Release
    add_test(
        NAME ogs_python_module
        COMMAND
            ${Python_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/ogs_python_module.py
            ${CMAKE_SOURCE_DIR}/Tests/Data/Parabolic/LiquidFlow/Flux/cube_1e3_calculatesurfaceflux.prj
    )
    set_tests_properties(
        ogs_python_module
        PROPERTIES LABELS "default" ENVIRONMENT_MODIFICATION
                   PYTHONPATH=path_list_append:$<TARGET_FILE_DIR:OpenGeoSys>
                   DISABLED $<NOT:$<TARGET_EXISTS:LiquidFlow>>
    )
endif()
# ---- Installation ----
install(TARGETS ogs RUNTIME DESTINATION bin)

set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "ogs"
                              "OGS Simulator"
)
cpack_add_component(
    ogs_cli
    DISPLAY_NAME "OGS THMC Simulator"
    DESCRIPTION "The command line interface for OpenGeoSys."
    GROUP Applications
)
