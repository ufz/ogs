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
    EXECUTABLES $<$<TARGET_EXISTS:ogs>:$<TARGET_FILE:ogs>> $<$<TARGET_EXISTS:DataExplorer>:$<TARGET_FILE:DataExplorer>>
    RESOLVED_DEPENDENCIES_VAR _r_deps
    UNRESOLVED_DEPENDENCIES_VAR _u_deps
    POST_EXCLUDE_REGEXES "/opt/local/lib/lib.*" # Disable macports zlib
  )
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
