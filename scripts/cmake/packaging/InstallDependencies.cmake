# Registers a target for installing its dependencies (dll / so files)
macro(InstallDependencies TARGET)
    set(INSTALL_DEPENDENCIES "${INSTALL_DEPENDENCIES};${TARGET}" CACHE INTERNAL "")
endmacro()
