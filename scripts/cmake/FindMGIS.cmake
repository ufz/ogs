# Find the MGIS includes and libraries
#
#    MGIS_INCLUDE_DIR - Where to find MGIS headers
#    MGIS_LIBRARY     - The MGIS library to link against.
#    MGIS_FOUND       - Do not attempt to use if "no" or undefined.
#
# MGIS, the MFront Generic Interface Support library
# See https://github.com/thelfer/MFrontGenericInterfaceSupport
#     http://tfel.sourceforge.net/generic-behaviours-interface.html

find_path(MGIS_INCLUDE_DIR MGIS)

find_library(MGIS_LIBRARY MFrontGenericInterface)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MGIS DEFAULT_MSG MGIS_LIBRARY MGIS_INCLUDE_DIR)
