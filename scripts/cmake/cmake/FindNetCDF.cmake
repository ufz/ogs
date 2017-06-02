# From https://github.com/jedbrown/cmake-modules
#
# Copyright $(git shortlog -s)
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice, this
#   list of conditions and the following disclaimer in the documentation and/or
#   other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
# - Find NetCDF
# Find the native NetCDF includes and library
#
#  NETCDF_INCLUDES    - where to find netcdf.h, etc
#  NETCDF_LIBRARIES   - Link these libraries when using NetCDF
#  NETCDF_FOUND       - True if NetCDF found including required interfaces (see below)
#
# Your package can require certain interfaces to be FOUND by setting these
#
#  NETCDF_CXX         - require the C++ interface and link the C++ library
#  NETCDF_F77         - require the F77 interface and link the fortran library
#  NETCDF_F90         - require the F90 interface and link the fortran library
#
# The following are not for general use and are included in
# NETCDF_LIBRARIES if the corresponding option above is set.
#
#  NETCDF_LIBRARIES_C    - Just the C interface
#  NETCDF_LIBRARIES_CXX  - C++ interface, if available
#  NETCDF_LIBRARIES_F77  - Fortran 77 interface, if available
#  NETCDF_LIBRARIES_F90  - Fortran 90 interface, if available
#
# Normal usage would be:
#  set (NETCDF_F90 "YES")
#  find_package (NetCDF REQUIRED)
#  target_link_libraries (uses_f90_interface ${NETCDF_LIBRARIES})
#  target_link_libraries (only_uses_c_interface ${NETCDF_LIBRARIES_C})

if (NETCDF_INCLUDES AND NETCDF_LIBRARIES)
  # Already in cache, be silent
  set (NETCDF_FIND_QUIETLY TRUE)
endif ()

find_path (NETCDF_INCLUDES netcdf.h
  HINTS NETCDF_DIR ENV NETCDF_DIR)

find_library (NETCDF_LIBRARIES_C       NAMES netcdf)
mark_as_advanced(NETCDF_LIBRARIES_C)

set (NetCDF_has_interfaces "YES") # will be set to NO if we're missing any interfaces
set (NetCDF_libs "${NETCDF_LIBRARIES_C}")

get_filename_component (NetCDF_lib_dirs "${NETCDF_LIBRARIES_C}" PATH)

macro (NetCDF_check_interface lang header libs)
  if (NETCDF_${lang})
    find_path (NETCDF_INCLUDES_${lang} NAMES ${header}
      HINTS "${NETCDF_INCLUDES}" NO_DEFAULT_PATH)
    find_library (NETCDF_LIBRARIES_${lang} NAMES ${libs}
      HINTS "${NetCDF_lib_dirs}" NO_DEFAULT_PATH)
    mark_as_advanced (NETCDF_INCLUDES_${lang} NETCDF_LIBRARIES_${lang})
    if (NETCDF_INCLUDES_${lang} AND NETCDF_LIBRARIES_${lang})
      list (INSERT NetCDF_libs 0 ${NETCDF_LIBRARIES_${lang}}) # prepend so that -lnetcdf is last
    else ()
      set (NetCDF_has_interfaces "NO")
      message (STATUS "Failed to find NetCDF interface for ${lang}")
    endif ()
  endif ()
endmacro ()

NetCDF_check_interface (CXX netcdfcpp.h netcdf_c++)
NetCDF_check_interface (F77 netcdf.inc  netcdff)
NetCDF_check_interface (F90 netcdf.mod  netcdff)

set (NETCDF_LIBRARIES "${NetCDF_libs}" CACHE STRING "All NetCDF libraries required for interface level")

# handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (NetCDF DEFAULT_MSG NETCDF_LIBRARIES NETCDF_INCLUDES NetCDF_has_interfaces)

mark_as_advanced (NETCDF_LIBRARIES NETCDF_INCLUDES)
