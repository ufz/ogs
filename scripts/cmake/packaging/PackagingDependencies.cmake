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

install(CODE "set(INSTALL_DIR \"${INSTALL_DIR}\")")
install(CODE "set(CMAKE_INSTALL_LIBDIR \"${CMAKE_INSTALL_LIBDIR}\")")
install(CODE "set(_rpath \"${_rpath}\")")
install(CODE "set(MKL_ROOT_DIR \"${MKL_ROOT_DIR}\")")
install(
    CODE [[
  file(GET_RUNTIME_DEPENDENCIES
    EXECUTABLES
        $<$<TARGET_EXISTS:ogs>:$<TARGET_FILE:ogs>>
        $<$<TARGET_EXISTS:DataExplorer>:$<TARGET_FILE:DataExplorer>>
        $<$<TARGET_EXISTS:testrunner>:$<TARGET_FILE:testrunner>>
        $<$<TARGET_EXISTS:RemoveGhostData>:$<TARGET_FILE:RemoveGhostData>>
    DIRECTORIES ${MKL_ROOT_DIR}/redist/intel64 ${MKL_ROOT_DIR}/../../tbb/latest/redist/intel64/vc_mt
    RESOLVED_DEPENDENCIES_VAR _r_deps
    UNRESOLVED_DEPENDENCIES_VAR _u_deps
    PRE_EXCLUDE_REGEXES "api-ms-" "ext-ms-"
    POST_EXCLUDE_REGEXES ".*system32/.*\\.dll"
  )
  find_program(PATCHELF_TOOL patchelf)
  foreach(_lib ${_r_deps})
    string(REGEX MATCH "libpetsc.*" _petsc_lib ${_lib})
    if(_petsc_lib AND EXISTS _ext/PETSc/lib/${_petsc_lib})
      if(PATCHELF_TOOL)
        execute_process(COMMAND patchelf --set-rpath ${_rpath} _ext/PETSc/lib/${_petsc_lib} COMMAND_ERROR_IS_FATAL ANY)
        message(STATUS "Patching RPATH of ${_petsc_lib} -> ${_rpath}")
      else()
          message(WARNING "patchelf tool not found: installed ogs binaries may not work (error: shared libraries not found)! "
            "Install the patchelf tool for proper runtime library search paths!")
      endif()
    endif()
  endforeach()
  foreach(_file ${_r_deps})
    file(INSTALL
      DESTINATION "${CMAKE_INSTALL_PREFIX}/${INSTALL_DIR}"
      TYPE SHARED_LIBRARY
      FILES "${_file}"
    )
  endforeach()
  list(LENGTH _u_deps _u_length)
  if("${_u_length}" GREATER 0)
    message(WARNING "Unresolved dependencies detected!\n${_u_deps}")
  endif()
]]
)
