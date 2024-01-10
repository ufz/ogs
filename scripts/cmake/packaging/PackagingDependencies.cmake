# Install dependencies via GET_RUNTIME_DEPENDENCIES.
if(NOT OGS_INSTALL_DEPENDENCIES)
    return()
endif()
if(WIN32)
    set(INSTALL_DIR ${CMAKE_INSTALL_BINDIR})
else()
    set(INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})
endif()
list(JOIN CMAKE_INSTALL_RPATH ":" _rpath)
install(CODE "set(CMAKE_BUILD_RPATH \"${CMAKE_BUILD_RPATH}\")")

install(CODE "set(INSTALL_DIR \"${INSTALL_DIR}\")")
install(CODE "set(CMAKE_INSTALL_LIBDIR \"${CMAKE_INSTALL_LIBDIR}\")")
install(CODE "set(_rpath \"${_rpath}\")")
file(TO_CMAKE_PATH "$ENV{MKLROOT}" MKL_ROOT_DIR)
install(CODE "set(MKL_ROOT_DIR \"${MKL_ROOT_DIR}\")")
install(
    CODE "set(OGS_INSTALL_DEPENDENCIES_PRE_EXCLUDES \"${OGS_INSTALL_DEPENDENCIES_PRE_EXCLUDES}\")"
)
install(
    CODE "set(OGS_INSTALL_DEPENDENCIES_POST_EXCLUDES \"${OGS_INSTALL_DEPENDENCIES_POST_EXCLUDES}\")"
)
install(
    CODE [[
  file(GET_RUNTIME_DEPENDENCIES
    EXECUTABLES
        $<$<TARGET_EXISTS:ogs>:$<TARGET_FILE:ogs>>
        $<$<TARGET_EXISTS:DataExplorer>:$<TARGET_FILE:DataExplorer>>
        $<$<TARGET_EXISTS:testrunner>:$<TARGET_FILE:testrunner>>
        $<$<TARGET_EXISTS:RemoveGhostData>:$<TARGET_FILE:RemoveGhostData>>
    DIRECTORIES
        ${MKL_ROOT_DIR}/redist/intel64
        ${MKL_ROOT_DIR}/../../tbb/latest/redist/intel64/vc_mt
        ${MKL_ROOT_DIR}/../../compiler/latest/windows/redist/intel64_win/compiler
        ${CMAKE_BUILD_RPATH}
    RESOLVED_DEPENDENCIES_VAR _r_deps
    UNRESOLVED_DEPENDENCIES_VAR _u_deps
    PRE_EXCLUDE_REGEXES "api-ms-" "ext-ms-" ${OGS_INSTALL_DEPENDENCIES_PRE_EXCLUDES}
    POST_EXCLUDE_REGEXES ".*system32/.*\\.dll" "/usr/lib/libtbb.so.12" ${OGS_INSTALL_DEPENDENCIES_POST_EXCLUDES}
  )
  find_program(PATCHELF_TOOL patchelf)
  foreach(_lib ${_r_deps})
    string(REGEX MATCH "libpetsc.*" _petsc_lib ${_lib})
    if(_petsc_lib AND EXISTS _ext/PETSc/lib/${_petsc_lib} AND NOT APPLE)
      if(PATCHELF_TOOL)
        message(STATUS "Patching RPATH of ${_petsc_lib} -> ${_rpath}")
        execute_process(COMMAND patchelf --set-rpath ${_rpath} _ext/PETSc/lib/${_petsc_lib} COMMAND_ERROR_IS_FATAL ANY)
      else()
          message(WARNING "patchelf tool not found: installed ogs binaries may not work (error: shared libraries not found)! "
            "Install the patchelf tool for proper runtime library search paths!")
      endif()
    endif()
  endforeach()
  message(STATUS "file(GET_RUNTIME_DEPENDENCIES) will install: ${_r_deps}")
  foreach(_file ${_r_deps})
    file(INSTALL
      DESTINATION "${CMAKE_INSTALL_PREFIX}/${INSTALL_DIR}"
      TYPE SHARED_LIBRARY
      FILES "${_file}"
      FOLLOW_SYMLINK_CHAIN
    )
  endforeach()
  list(LENGTH _u_deps _u_length)
  if("${_u_length}" GREATER 0)
    message(WARNING "Unresolved dependencies detected!\n${_u_deps}")
  endif()
]]
)
