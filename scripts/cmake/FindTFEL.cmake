# Usage:
#
# # TFELHOME optional for Conan or custom install locations
# set(TFELHOME ${CONAN_TFEL_ROOT})
# find_package(TFEL REQUIRED)
#
# Once found, declares these variables:
#   - TFEL_FOUND
#   - TFEL_INCLUDE_PATH
#   - TFEL_LIBRARIES
#   - TFEL_CXX_STANDARD
#   - TFEL_VERSION

#if(DEFINED ENV{TFELHOME})
#    set(TFELHOME $ENV{TFELHOME})
#endif()
message(STATUS "tfelhome: ${TFELHOME}")

foreach(tool mfront tfel-check tfel-config mfront-query)
    string(TOUPPER ${tool} toolVar)
    string(REPLACE "-" "_" toolVar ${toolVar})
    find_program(${toolVar} ${tool} HINTS ${TFELHOME} PATH_SUFFIXES bin)
endforeach()

IF(TFEL_CONFIG AND MFRONT AND MFRONT_QUERY)

  execute_process(COMMAND ${TFEL_CONFIG} "--cxx-standard"
      RESULT_VARIABLE TFEL_CXX_STANDARD_AVAILABLE
      OUTPUT_VARIABLE TFEL_CXX_STANDARD
      OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT TFEL_CXX_STANDARD_AVAILABLE EQUAL 0)
      set(TFEL_CXX_STANDARD 11)
  endif()

  execute_process(COMMAND ${TFEL_CONFIG} "--version"
      OUTPUT_VARIABLE TFEL_VERSION_FULL
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
  message(STATUS "TFEL_VERSION_FULL: ${TFEL_VERSION_FULL}")
  string(REPLACE " " ";" TFEL_VERSION ${TFEL_VERSION_FULL})
  list(GET TFEL_VERSION 1 TFEL_VERSION)

  set(TFEL_PYTHON_BINDINGS OFF)
  execute_process(COMMAND ${TFEL_CONFIG} "--python-version"
      RESULT_VARIABLE TFEL_PYTHON_BINDINGS
      OUTPUT_VARIABLE TFEL_PYTHON_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
  )
  if(TFEL_PYTHON_BINDINGS EQUAL 0)
      set(TFEL_PYTHON_BINDINGS ON)
  endif()

  find_path(TFEL_CONFIG_INCLUDE_PATH TFELConfig.hxx
      HINTS ${TFELHOME}/include/TFEL/Config
  )
  get_filename_component(TFEL_INCLUDE_PATH
      ${TFEL_CONFIG_INCLUDE_PATH}/../.. ABSOLUTE CACHE
  )

  set(tfel_libs
      TFELTests
      TFELException
      TFELUtilities
      TFELMaterial
      TFELMath
      MTestFileGenerator
  )
  if(${TFEL_CXX_STANDARD} LESS 17)
      list(APPEND tfel_libs TFELPhysicalConstants)
  endif()

  foreach(lib ${tfel_libs})
      find_library(${lib}_LIBRARY ${lib} HINTS ${TFELHOME} PATH_SUFFIXES lib)
      list(APPEND TFEL_LIBRARIES ${${lib}_LIBRARY})
  endforeach(lib)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(TFEL DEFAULT_MSG
      TFEL_VERSION
      TFEL_INCLUDE_PATH
      TFEL_LIBRARIES
      TFEL_CXX_STANDARD
  )

  set(TFEL_FOUND ON)

else(TFEL_CONFIG AND MFRONT AND MFRONT_QUERY)

  set(TFEL_FOUND OFF)

endif(TFEL_CONFIG AND MFRONT AND MFRONT_QUERY)
