# Install dependencies via GET_RUNTIME_DEPENDENCIES.
if(NOT OGS_INSTALL_DEPENDENCIES)
    return()
endif()
install(
    CODE [[
  include(GNUInstallDirs)
  if(WIN32)
    set(INSTALL_DIR ${CMAKE_INSTALL_FULL_BINDIR})
  else()
    set(INSTALL_DIR ${CMAKE_INSTALL_FULL_LIBDIR})
  endif()
  file(GET_RUNTIME_DEPENDENCIES
    EXECUTABLES
        $<$<TARGET_EXISTS:ogs>:$<TARGET_FILE:ogs>>
        $<$<TARGET_EXISTS:DataExplorer>:$<TARGET_FILE:DataExplorer>>
        $<$<TARGET_EXISTS:testrunner>:$<TARGET_FILE:testrunner>>
    RESOLVED_DEPENDENCIES_VAR _r_deps
    UNRESOLVED_DEPENDENCIES_VAR _u_deps
    POST_EXCLUDE_REGEXES "/opt/local/lib/lib.*" # Disable macports zlib
  )
  find_program(PATCHELF_TOOL patchelf)
  if(PATCHELF_TOOL)
    foreach(_lib ${_r_deps})
      string(REGEX MATCH "libpetsc.*" _petsc_lib ${_lib})
      if(_petsc_lib AND EXISTS _ext/PETSc/lib/${_petsc_lib})
        message(STATUS "Patching RPATH of ${_petsc_lib}")
        execute_process(COMMAND patchelf --set-rpath '$ORIGIN:$ORIGIN/../lib:/usr/local/openmpi/lib:/usr/lib/openmpi' _ext/PETSc/lib/${_petsc_lib})
      endif()
    endforeach()
  endif()
  file(INSTALL ${_r_deps}
    DESTINATION ${INSTALL_DIR}
    FOLLOW_SYMLINK_CHAIN
  )
  list(LENGTH _u_deps _u_length)
  if("${_u_length}" GREATER 0)
    message(WARNING "Unresolved dependencies detected!\n${_u_deps}")
  endif()
]]
)
